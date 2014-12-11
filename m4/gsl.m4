# AC_SEARCH_GSL(actionIfFound, actionIfNotFound)
AC_DEFUN([AC_SEARCH_GSL], [

AC_PATH_PROG(GSLCONFIG, gsl-config, [], [$PATH])
if test -f "$GSLCONFIG"; then
  GSL_CPPFLAGS=`$GSLCONFIG --cflags`
  GSL_LDFLAGS=`$GSLCONFIG --libs`
else
  AC_MSG_ERROR([GSL cannot be found!])
  exit 1
fi
AC_SUBST(GSL_CPPFLAGS)
AC_SUBST(GSL_LDFLAGS)
$1
])

