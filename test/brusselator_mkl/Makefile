# These override the defaults in ridc.mk
SRCDIR  = ../../src
OBJDIR  = ../../build

# Path to final executable
EXE     = brusselator.exe

# Allows us to just run 'make' (this needs to be above the include)
default: brusselator_mkl

# Include the build rules from ridc.mk
include $(SRCDIR)/ridc.mk

.PHONY: brusselator_mkl

# Add this directory to include path (for ode.h) and MKL libs
brusselator_mkl: CPPFLAGS+=-DIMPLICIT -I$(realpath .) -lmkl_intel_lp64 -lmkl_core -lmkl_sequential -lpthread -lm

# Build the executable with $(OBJECTS) defined in ridc.mk
brusselator_mkl: $(OBJECTS)
	$(Q)$(CXX) -o $(EXE) $^  brusselator.cpp $(CPPFLAGS)
	$(Q)$(SAY) "CXX    $(EXE)"

