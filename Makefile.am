ACLOCAL_AMFLAGS = -I m4
AUTOMAKE_OPTIONS = foreign
SUBDIRS = src . \
	examples/explicit \
	examples/implicit \
	examples/implicit_mkl \
	examples/brusselator_mkl\
	examples/brusselator_gsl

# DOXYGEN SUPPORT
include $(top_srcdir)/aminclude.am

# ensure the distribution of the doxygen configuration file
EXTRA_DIST = Doxyfile.in
