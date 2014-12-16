AC_PREREQ([2.69])
AC_INIT([compressor],[1.0.0],[stefano.carrazza@mi.infn.it])
AC_CONFIG_SRCDIR([src/main.cc])
AC_CONFIG_HEADERS([inc/config.h])
AC_CONFIG_MACRO_DIR([m4])
AM_INIT_AUTOMAKE([subdir-objects])
m4_ifdef([AM_SILENT_RULES], [AM_SILENT_RULES([yes])])
AM_MAINTAINER_MODE([disable])

# checks for programs
AC_PROG_CXX
AC_PROG_INSTALL

# check for root
#AC_SEARCH_ROOT
AC_SEARCH_LHAPDF
AC_SEARCH_GSL

AC_CONFIG_FILES([Makefile inc/Makefile src/Makefile])
AC_OUTPUT
