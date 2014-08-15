OBJDIR  ?= ../build
SRCDIR  ?= ../src
TESTDIR ?= ../test

sources = ridc.cpp driver.cpp ridc_utils.cpp
SOURCES = $(addprefix $(SRCDIR)/, $(sources))
OBJECTS = $(SOURCES:$(SRCDIR)/%.cpp=$(OBJDIR)/%.o)

EXE ?= $(BINDIR)/ridc.exe

# Compiler flags
CPPFLAGS += -fopenmp -g

# determine vebosity level
ifdef V
	ifeq ("$(origin V)", "command line")
		VERBOSE=$(V)
	endif
else
	VERBOSE=0
endif

# Set command prefixes based on verbosity
ifeq ($(VERBOSE), 1)
	SAY=@true
	Q=
else
	SAY=echo
	Q=@
endif

# if make -s (silent mode) suppress all output
ifneq ($(findstring s, $(MAKEFLAGS)),)
	Q=@
	SAY=@true
endif

$(OBJECTS): $(OBJDIR)/%.o: $(SRCDIR)/%.cpp $(OBJDIR) $(BINDIR)
	$(Q)$(CXX) -c -o $@ $< $(CPPFLAGS)
	$(Q)$(SAY) "CXX    $@"
$(OBJDIR):
	$(Q)mkdir -p $@

.PHONY: clean
clean:
	$(Q)$(RM) $(OBJDIR)/*.o $(EXE)
	$(Q)$(SAY) "Build directory clean"
