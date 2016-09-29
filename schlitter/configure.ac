#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.
AC_INIT(schlitter, 1.0.0, jens@jkleinj.eu)

AC_PREREQ([2.63])
AC_CONFIG_SRCDIR([src/arg.c])
AM_INIT_AUTOMAKE([foreign -Wall -Werror])

# Checks for programs.
AC_PROG_CC_C_O
AC_PROG_INSTALL

# Checks for libraries.
AC_CHECK_LIB([m], [main], [],
				AC_MSG_ERROR([Could not find required C math library.]))

# Checks for header files.
AC_CHECK_HEADERS([float.h stdlib.h string.h])

# Checks for typedefs, structures, and compiler characteristics.
AC_TYPE_SIZE_T

# Checks for library functions.
AC_FUNC_MALLOC
AC_FUNC_REALLOC
AC_CHECK_FUNCS([pow sqrt])

AC_CHECK_LIB([gslcblas], main)
AC_CHECK_LIB([gsl], main)

# GSL
AC_ARG_WITH(gsl, [  --with-gsl=<path>      prefix of GSL installation],
                 [  CPPFLAGS="$CPPFLAGS -I $withval/include"
                    LDFLAGS="$LDFLAGS -L $withval/lib"])

AC_CHECK_HEADER([gsl/gsl_math.h], [],
             AC_MSG_ERROR([Could not find required GSL header files.
                           Consider using --with-gsl=<my-path-to-gsl>
                           or installing the GSL if required.
                           Read the INSTALL file for details.]))

AC_CHECK_LIB([gslcblas], [main], [],
            AC_MSG_ERROR([Could not find required GSL BLAS library.
                          Consider using --with-gsl=<my-path-to-gsl>
                          or installing the GSL if required.
                          Read the INSTALL file for details.]))

AC_CHECK_LIB([gsl], [main], [],
             AC_MSG_ERROR([Could not find required GSL library.
                          Consider using --with-gsl=<my-path-to-gsl>
                          or installing the GSL if required.
                          Read the INSTALL file for details.]))

# enable debug
AC_ARG_ENABLE(debug, [  --enable-debug  Turn on debugging], [ac_use_debug="yes"], [ac_use_debug="no"])
if test "$ac_use_debug" = "yes"; then
    CFLAGS="-ggdb"
    AC_DEFINE(DEBUG,1,[Define to compile with DEBUG settings.])
fi

# enable profiling 
AC_ARG_ENABLE(profiling, [  --enable-profiling  Turn on profiling], [ac_use_profiling="yes"], [ac_use_profiling="no"])
if test "$ac_use_profiling" = "yes"; then
    CFLAGS="-pg -O2"
    AC_DEFINE(PROFILING,1,[Define to compile with PROFILING settings.])
fi

# Output files
AC_CONFIG_FILES([Makefile src/Makefile tests/Makefile])
AC_CONFIG_HEADERS([src/config.h])
AC_OUTPUT

