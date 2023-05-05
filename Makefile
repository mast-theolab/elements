FC := gfortran

SRCDIR = src
BUILDIR = build

vpath %.f90 src/lib:src/prog

# List 
# The order is IMPORTANT! To handle inter-dependencies
SRC_LIB = atomic.f90 exception.f90 math.f90 physic.f90 string.f90 output.f90 \
	parsefchk.f90 basisset.f90 electronic.f90 moldata.f90 transdata.f90 \
	input.f90 exc_sos.f90

SRC_GMCD = gmcd_output.f90 gmcd_legacy.f90 mcd_tensor.f90

OBJ_LIB = $(patsubst %.f90, $(BUILDIR)/%.o, $(SRC_LIB))
OBJ_GMCD = $(patsubst %.f90, $(BUILDIR)/%.o, $(SRC_GMCD))

EXE_GMCD = $(basename $(lastword $(SRC_GMCD)))

.PHONY: clean docs

# Compilation mode
MODE := run
ifeq ($(MODE), debug)
	LDFLAGS = -g
	FCFLAGS = -g -c -Wall -Wextra -Wconversion -Og -pedantic -fcheck=bounds -fmax-errors=5
else
	FCFLAGS = -O4 -c
	LDFLAGS = 
endif
FCFLAGS += -J$(BUILDIR)

# This must be the first rule
all: $(EXE_GMCD)

$(EXE_GMCD): $(OBJ_LIB) $(OBJ_GMCD)
	cd $(BUILDIR) && \
	$(FC) $(LDFLAGS) -o $@ $(patsubst $(BUILDIR)/%, %, $^) && \
	cp $@ ../bin

clean:
	-rm -rf build/

docs:
	ford elements.md

# Compiler steps for all objects
# The pipe ensures that the directory is present
$(OBJ_LIB) : $(BUILDIR)/%.o : %.f90 | $(BUILDIR)
	$(FC) $(FCFLAGS) -o $@ $<

$(OBJ_GMCD) : $(BUILDIR)/%.o : %.f90 | $(BUILDIR)
	$(FC) $(FCFLAGS) -o $@ $<

$(BUILDIR):
	mkdir -p $(BUILDIR)

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
