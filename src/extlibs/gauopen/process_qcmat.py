"""A simple script to parse and process GauOpen

The GauOpen library provides 2 (or 3) low-level Fortran files:
* qcmatrix.F
* qcmatrixio.F
* qcmatrixiou.F since version 3
that can produce two sets of routines with the same names but different
kinds for the integers.
This is not a problem if separate executables are generated but becomes
complicated if we want to integrate all in a single space.
To avoid this, we create two duplicates of each file, one based on
32-bits integers (_32) and one based on 64-bits (_64).

For simplicity, we work in two passes:

1. build a database of all defined routines over both libraries and
   act as a very simple preprocessor for USE_I8 (extensible)
2. fix any reference to the routines in calls, to ensure there is
   no overlap.

Since the files are not expected to be large, we work in memory for
simplicity.
We keep a memory check to ensure we do not go out of memory.
"""

import os
import sys
import re


SRC_FILES_BASE = ("qcmatrix.F", "qcmatrixio.F")
SRC_FILES_EXTRA = ("qcmatrixiou.F", )
MAX_MEM = 1_000_000


def split_line(line: str) -> list[str]:
    """Split line based on basic parsing conventions."""
    newline = []
    if line[0].lower() in ('*', 'c'):
        if len(line) > 80:
            if ' ' not in line:
                print('ERROR: Could not find one space in comment line:')
                print(f'       {line}')
                print('       Something is strange, stopping.')
                sys.exit(2)
            line1, line2 = line.rsplit(' ', maxsplit=1)
            i = 1
            while line1[i] == ' ':
                i += 1
            prefix = 'C' + (i-1)*' '
            newline.append(f'{line1}\n')
            newline.append(f'{prefix}{line2}')
        else:
            newline.append(line)
    else:
        for char in (',', '(', ' '):
            if 0 <= line[:72].rfind(char) <= 72:
                i = line[:72].rfind(char)
                line1 = line[:i+1].rstrip()
                line2 = line[i+1:]
                if line1[5] != ' ':
                    i = 6
                    while line1[i] == ' ':
                        i += 1
                    prefix = 5*' ' + line1[5] + (i-6)*' '
                else:
                    prefix = 5*' ' + '$' + 2*' '
                newline.append(f'{line1}\n')
                newline.append(f'{prefix}{line2}')
                break
        else:
            print('ERROR: Could not find how to properly cut line:')
            print(f'       {newline}')
            print('      Cowardly quitting.')
            sys.exit(2)

    return newline


def main():
    """Run the main script."""
    # First, let us check that all base files exist and check for extras
    # Check memory consumption as well.
    src_files = []
    estimated_mem = 0
    for file in SRC_FILES_BASE:
        if not os.path.exists(file):
            print(f"ERROR: File {file} not found.  Cannot proceed.")
        src_files.append(file)
        estimated_mem += os.path.getsize(file)
    for file in SRC_FILES_EXTRA:
        if os.path.exists(file):
            src_files.append(file)
            estimated_mem += os.path.getsize(file)

    # We plan 3 times the memory, the source, 32 and 64 bits versions.
    print(f"Total expected memory (in kB): {3*estimated_mem/1000:.0f}")
    if estimated_mem > MAX_MEM:
        print("ERROR: Memory requirement too large.")
        print("       Increase MAX_MEM if this is still acceptable.")
        sys.exit(1)

    lines_src = {}
    for file in src_files:
        with open(file, "r", encoding="utf-8") as fobj:
            lines_src[file] = fobj.readlines()

    # Extract list of routines
    # We take advantage of Gaussian standard structure:
    # *Deck <Routine>
    # <Declaration>
    # We assume that there is no name overlap between variables and routines
    routines = []
    for file, lines in lines_src.items():
        for line in lines:
            if line.lower().startswith('*deck'):
                routines.append(line.split()[-1])
    if 'Close_MatF' in routines:
        gauopen_version = 2
    else:
        gauopen_version = 3
    # Now let us build a general regexp to find all cases
    key_routines = re.compile(f"\\b({'|'.join(routines)})\\b", flags=re.I)
    key_cppI64 = '#ifdef USE_I8'

    # Now build the lines
    lines_i32 = {}
    lines_i64 = {}
    for file, lines in lines_src.items():
        lines_i32[file] = []
        lines_i64[file] = []
        in_define = False
        for line in lines:
            if res := key_routines.findall(line):
                if len(res) == 1:
                    newline32 = line.replace(res[0], res[0]+'_I4')
                    newline64 = line.replace(res[0], res[0]+'_I8')
                else:
                    newline32 = key_routines.sub(
                        lambda obj: obj.group(0) + '_I4', line)
                    newline64 = key_routines.sub(
                        lambda obj: obj.group(0) + '_I8', line)
                # Check if substitution led to excessive length for fixed
                # Fortran.
                # Since this would be true for both 32 and 64 bits versions,
                # do the check on one.
                if len(newline32.rstrip('\n')) > 72:
                    lines_i32[file].extend(split_line(newline32))
                    lines_i64[file].extend(split_line(newline64))
                else:
                    lines_i32[file].append(newline32)
                    lines_i64[file].append(newline64)
            elif line.strip() == key_cppI64:
                in_define = 64
            elif in_define:
                match (line.strip()):
                    case '#else':
                        in_define = 32
                    case '#endif':
                        in_define = False
                    case _:
                        if in_define == 64:
                            lines_i64[file].append(line)
                        else:
                            lines_i32[file].append(line)
            else:
                lines_i32[file].append(line)
                lines_i64[file].append(line)

    # Now let us build the files
    for file in lines_src:
        base = os.path.splitext(file)[0]
        with open(f"{base}4.F", "w", encoding="utf-8") as fobj32, \
                open(f"{base}8.F", "w", encoding="utf-8") as fobj64:
            fobj32.writelines(lines_i32[file])
            fobj64.writelines(lines_i64[file])

    with open('gauopen.version', 'w', encoding='utf-8') as fobj:
        fobj.write(f'{gauopen_version}\n')


if __name__ == '__main__':
    main()
