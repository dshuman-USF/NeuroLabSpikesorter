#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

AC_PREREQ(2.57)
AC_INIT(USF Spikesorter, 0.18.3, dshuman@usf.edu, spikesorter)
AM_INIT_AUTOMAKE
AC_CONFIG_SRCDIR([cdf.c])
AC_CONFIG_HEADER([config.h])
AC_CANONICAL_HOST

# Checks for programs.
AC_PROG_CXX
AC_PROG_CC
AC_PROG_LN_S
AC_C_BIGENDIAN

# Checks for libraries.
AC_SEARCH_LIBS(sem_post, pthread)
AC_SEARCH_LIBS(floor, m)
AC_SEARCH_LIBS(ao_play, ao)
PKG_CHECK_MODULES([FFTW3], [fftw3])
PKG_CHECK_MODULES([GSL], [gsl])
PKG_CHECK_MODULES([X], [x11])
PKG_CHECK_MODULES([GTK], [gtk+-3.0])

# Checks for header files.
AC_HEADER_STDC
AC_CHECK_HEADERS([limits.h malloc.h stdlib.h string.h unistd.h])

# Checks for typedefs, structures, and compiler characteristics.
AC_HEADER_STDBOOL
AC_C_CONST
AC_C_INLINE
AC_TYPE_SIZE_T

# Checks for library functions.
AC_FUNC_MEMCMP
AC_FUNC_REALLOC
AC_CHECK_FUNCS([floor memmove memset pow sqrt strdup strrchr ao_play])
AC_REPLACE_FUNCS(asprintf)
AC_CONFIG_FILES([Makefile])
AC_CONFIG_FILES([spikesort_control_panel],[chmod 755 spikesort_control_panel])
AC_CONFIG_FILES([alert],[chmod 755 alert])
AC_CONFIG_FILES([tkss.tcl],[chmod 755 tkss.tcl])
AC_OUTPUT
