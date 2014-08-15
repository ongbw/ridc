# These override the defaults in ridc.mk
SRCDIR  = ../../src
OBJDIR  = ../../build

# Path to final executable
EXE     = explicit_ode1.exe

# Allows us to just run 'make' (this needs to be above the include)
default: explicit_ode1

# Include the build rules from ridc.mk
include $(SRCDIR)/ridc.mk

.PHONY: explicit_ode1

# Add this directory to include path (for ode.h)
explicit_ode1: CPPFLAGS+=-DEXPLICIT -I$(realpath .)

# Build the executable with $(OBJECTS) defined in ridc.mk
explicit_ode1: $(OBJECTS)
	$(Q)$(CXX) -o $(EXE) $^ ode.cpp $(CPPFLAGS)
	$(Q)$(SAY) "CXX    $(EXE)"

