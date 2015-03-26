# Process this file with autoconf to produce a configure script.

AC_INIT([libridc], [0.1], [ongbw@mtu.edu])
AC_CONFIG_SRCDIR([src/ridc.cpp])
AM_INIT_AUTOMAKE([1.9])
AC_PROG_LIBTOOL
LT_INIT

# use the C++ compiler for the following checks
AC_LANG([C++])

# Checks for programs.
AC_PROG_CXX
AC_OPENMP

# distribute additional compiler and linker flags
# --> set these variables instead of CXXFLAGS or LDFLAGS
AC_SUBST([AM_CXXFLAGS])
AC_SUBST([AM_LDFLAGS])
AC_SUBST([LIBS])

# files to generate via autotools (.am or .in source files)
AC_CONFIG_FILES([Makefile])
AC_CONFIG_FILES([src/Makefile])
AC_CONFIG_FILES([examples/explicit/Makefile])
AC_CONFIG_FILES([examples/implicit/Makefile])
AC_CONFIG_FILES([examples/implicit_mkl/Makefile])
AC_CONFIG_FILES([examples/brusselator_mkl/Makefile])
AC_CONFIG_FILES([examples/brusselator_gsl/Makefile])


######################################################################
# DOXYGEN SUPPORT
######################################################################

DX_HTML_FEATURE(ON)
DX_CHM_FEATURE(OFF)
DX_CHI_FEATURE(OFF)
DX_MAN_FEATURE(OFF)
DX_RTF_FEATURE(OFF)
DX_XML_FEATURE(OFF)
DX_PDF_FEATURE(OFF)
DX_PS_FEATURE(OFF)

AC_CONFIG_FILES([Doxyfile])
DX_INIT_DOXYGEN([$PACKAGE_NAME],[Doxyfile])

######################################################################

AC_OUTPUT
