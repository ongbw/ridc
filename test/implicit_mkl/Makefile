# These override the defaults in ridc.mk
SRCDIR  = ../../src
OBJDIR  = ../../build
TESTDIR = ../../test

# Path to final executable
EXE     = implicit_ode_mkl.exe

# Allows us to just run 'make' (this needs to be above the include)
default: implicit_ode_mkl

# Include the build rules from ridc.mk
include $(SRCDIR)/ridc.mk

.PHONY: implicit_ode_mkl

# Add this directory to include path (for ode.h)
implicit_ode_mkl: CPPFLAGS+=-DIMPLICIT -I$(realpath .) -lmkl_intel_lp64 -lmkl_core -lmkl_sequential -lpthread -lm


# Build the executable with $(OBJECTS) defined in ridc.mk
implicit_ode_mkl: $(OBJECTS)
	$(Q)$(CXX) -o $(EXE) $^ ode_implicit.cpp $(CPPFLAGS)
	$(Q)$(SAY) "CXX    $(EXE)"

