FC := gfortran

SRCDIR = src
BUILDIR = build
BINDIR = bin

vpath %.f90 src/lib:src/prog

# List 
# The order is IMPORTANT! To handle inter-dependencies
SRC_LIB_CMDARG = parse_cmdline.f90 parse_cmdline_parser.f90 \
	parse_cmdline_addarg.f90 parse_cmdline_setval.f90 parse_cmdline_getval.f90
SRC_LIB = numeric.f90 atomic.f90 exception.f90 math.f90 physic.f90 string.f90 \
	output.f90 $(SRC_LIB_CMDARG) parsefchk.f90 basisset.f90 electronic.f90 \
	moldata.f90 transdata.f90 input.f90 exc_sos.f90

SRC_GMCD = gmcd_output.f90 gmcd_legacy.f90 mcd_tensor.f90
SRC_DEVARG = test_cmdline.f90

OBJ_LIB = $(patsubst %.f90, $(BUILDIR)/%.o, $(SRC_LIB))
OBJ_GMCD = $(patsubst %.f90, $(BUILDIR)/%.o, $(SRC_GMCD))
OBJ_DEVARG = $(patsubst %.f90, $(BUILDIR)/%.o, $(SRC_DEVARG))
OBJ_PROGS = $(OBJ_GMCD) $(OBJ_DEVARG)

EXE_GMCD = $(basename $(lastword $(SRC_GMCD)))
EXE_DEVARG = $(basename $(lastword $(SRC_DEVARG)))

.PHONY: clean docs

# Compilation mode
MODE := run
ifeq ($(MODE), debug)
	LDFLAGS = -g
	FCFLAGS = -g -c -Wall -Wextra -Wconversion -Og -pedantic -fcheck=bounds -fmax-errors=5 -std=f2008
else
	FCFLAGS = -O4 -c
	LDFLAGS = 
endif
FCFLAGS += -J$(BUILDIR)

# This must be the first rule
all: $(EXE_GMCD) | $(BINDIR)

$(EXE_GMCD): $(OBJ_LIB) $(OBJ_GMCD)
	cd $(BUILDIR) && \
	$(FC) $(LDFLAGS) -o $@ $(patsubst $(BUILDIR)/%, %, $^) && \
	cp $@ ../$(BINDIR)

clean:
	-rm -rf $(BUILDIR)

docs:
	ford elements.md

# Compiler steps for all objects
# The pipe ensures that the directory is present
$(OBJ_LIB) : $(BUILDIR)/%.o : %.f90 | $(BUILDIR)
	$(FC) $(FCFLAGS) -o $@ $<

$(OBJ_PROGS) : $(BUILDIR)/%.o : %.f90 | $(BUILDIR)
	$(FC) $(FCFLAGS) -o $@ $<

$(BUILDIR):
	mkdir -p $(BUILDIR)

$(BINDIR):
	mkdir -p $(BINDIR)

# # Linker
# $(PRG_NAME): $(OBJS)
# 	cd $(BUILDIR) ; \
# 	$(FC) $(LDFLAGS) -o $@ $(patsubst $(BUILDIR)/%,%,$^)

# # Rule for main program, which depends on all modules
# $(PRG_OBJ): $(MOD_OBJS)

debug:
	@echo "SRCS = $(SRCS)"
	@echo "OBJS = $(OBJS)"
	@echo "MODS = $(MODS)"
	@echo "MOD_OBJS = $(MOD_OBJS)"
	@echo "PROGRAM = $(PRG_NAME)"
	@echo "PRG_OBJ = $(PRG_OBJ)"

tests: test_args

test_args: $(EXE_DEVARG)

$(EXE_DEVARG): $(OBJ_LIB) $(OBJ_DEVARG) | $(BINDIR)
	cd $(BUILDIR) && \
	$(FC) $(LDFLAGS) -o $@ $(patsubst $(BUILDIR)/%, %, $^) && \
	cp $@ ../$(BINDIR)
