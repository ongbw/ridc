# These override the defaults in ridc.mk
SRCDIR  = ../../src
OBJDIR  = ../../build

# Path to final executable
EXE     = implicit_ode1.exe

# Allows us to just run 'make' (this needs to be above the include)
default: implicit_ode1

# Include the build rules from ridc.mk
include $(SRCDIR)/ridc.mk

.PHONY: implicit_ode1

# Add this directory to include path (for ode.h)
implicit_ode1: CPPFLAGS+=-DIMPLICIT -I$(realpath .)

# Build the executable with $(OBJECTS) defined in ridc.mk
implicit_ode1: $(OBJECTS)
	$(Q)$(CXX) -o $(EXE) $^ ode_implicit.cpp $(CPPFLAGS)
	$(Q)$(SAY) "CXX    $(EXE)"

