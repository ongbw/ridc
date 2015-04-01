LIBS = -lmkl_intel_lp64 -lmkl_sequential -lmkl_core -lm -lpthread
OBJS = \
	brusselator.o \
	radau.o

all: brusselator_radau.exe

brusselator.o: brusselator.cpp

radau.o : radau.cpp

brusselator_radau.exe: $(OBJS)
	$(CXX) $(CPPFLAGS) $(OBJS) -o $@ $(LIBS)
