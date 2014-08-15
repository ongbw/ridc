# These override the defaults in ridc.mk
SRCDIR  = ../../src
OBJDIR  = ../../build

# Path to final executable
EXE     = brusselator.exe

# Allows us to just run 'make' (this needs to be above the include)
default: brusselator_gsl

# Include the build rules from ridc.mk
include $(SRCDIR)/ridc.mk

.PHONY: brusselator_gsl

# Add this directory to include path (for ode.h), appropriate libraries
# NOTE: some GSL installations will require an explicit link to -lgslcblas:
# brusselator_gsl: CPPFLAGS+=-DIMPLICIT -I$(realpath .) -lgsl -lgslcblas -lm 

brusselator_gsl: CPPFLAGS+=-DIMPLICIT -I$(realpath .) -lgsl -lm

# Build the executable with $(OBJECTS) defined in ridc.mk
brusselator_gsl: $(OBJECTS)
	$(Q)$(CXX) -o $(EXE) $^  brusselator.cpp $(CPPFLAGS)
	$(Q)$(SAY) "CXX    $(EXE)"

