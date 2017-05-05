AC_INIT([parallel], [20170422], [bug-parallel@gnu.org])
AM_INIT_AUTOMAKE([-Wall -Werror foreign])
AC_CONFIG_HEADERS([config.h])
AC_CONFIG_FILES([
  Makefile
  src/Makefile
])

AC_ARG_ENABLE(documentation,
  AC_HELP_STRING([--disable-documentation],
		 [Omit building and installing the documentation. (default=no)]),,
		 [enable_documentation=yes])
AM_CONDITIONAL([DOCUMENTATION], [test x$enable_documentation = xyes])
AC_PROG_LN_S
AC_OUTPUT
