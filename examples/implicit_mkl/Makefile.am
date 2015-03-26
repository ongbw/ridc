AM_CXXFLAGS = --pedantic -Wall $(OPENMP_CXXFLAGS)
AM_CPPFLAGS = $(CPPFLAGS) -I$(top_srcdir)/src 


LOG_DRIVER = $(top_srcdir)/examples/implicit_mkl/run.sh

TESTS = implicit_mkl
check_PROGRAMS = implicit_mkl
implicit_mkl_SOURCES = implicit.cpp  driver.cpp
implicit_mkl_LDADD = $(top_builddir)/src/libridc.la -lmkl_intel_lp64 -lmkl_core -lmkl_sequential -lpthread -lm
