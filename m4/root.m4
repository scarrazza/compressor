# AC_SEARCH_ROOT(actionIfFound, actionIfNotFound)
AC_DEFUN([AC_SEARCH_ROOT], [

AC_PATH_PROG(ROOTCONFIG, root-config, [], [$PATH])
if test -f "$ROOTCONFIG"; then
  ROOT_CPPFLAGS=`$ROOTCONFIG --cflags`
  ROOT_LDFLAGS=`$ROOTCONFIG --glibs`
else
  AC_MSG_ERROR([ROOT cannot be found!])
  exit 1
fi
AC_SUBST(ROOT_CPPFLAGS)
AC_SUBST(ROOT_LDFLAGS)
$1
])

