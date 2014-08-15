# These override the defaults in ridc.mk
SRCDIR  = ../../src
OBJDIR  = ../../build

# Path to final executable
EXE     = reaction_diffusion.exe

# Allows us to just run 'make' (this needs to be above the include)
default: reaction_diffusion

# Include the build rules from ridc.mk
include $(SRCDIR)/ridc.mk

.PHONY: reaction_diffusion

# Add this directory to include path (for ode.h)
reaction_diffusion: CPPFLAGS+=-DEXPLICIT -I$(realpath .)

# Build the executable with $(OBJECTS) defined in ridc.mk
reaction_diffusion: $(OBJECTS)
	$(Q)$(CXX) -o $(EXE) $^ reaction_diffusion.cpp $(CPPFLAGS)
	$(Q)$(SAY) "CXX    $(EXE)"

