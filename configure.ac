# -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

AC_PREREQ(2.57)
AC_INIT([StatPatternRecognition],[3.3.2],[andy.buckley@durham.ac.uk],[SPR])
AC_CONFIG_SRCDIR([src/SprRandomNumber.cc])
#AC_CONFIG_HEADER([include/StatPatternRecognition.hh])

AM_INIT_AUTOMAKE(dist-bzip2)
AC_LANG(C++)

# Checks for programs.
AC_PROG_CXX
AC_PROG_CC
AC_PROG_CPP
AC_PROG_INSTALL
AC_PROG_LN_S
AC_PROG_MAKE_SET
AC_PROG_LIBTOOL

# Checks for header files.
AC_HEADER_STDC
AC_CHECK_HEADERS([float.h stdlib.h sys/time.h unistd.h])

# Checks for typedefs, structures, and compiler characteristics.
AC_HEADER_STDBOOL
AC_CHECK_HEADERS([stddef.h])
AC_C_CONST
AC_C_INLINE
AC_HEADER_TIME
AC_TYPE_SIZE_T

## Checks for library functions.
AC_FUNC_ERROR_AT_LINE
AC_FUNC_STAT
AC_CHECK_FUNCS([gettimeofday pow sqrt floor])

## Default compiler flags
DEFAULTCXXFLAGS="-O4"

## ROOT installation
AC_ARG_VAR([ROOTPATH], [path to ROOT installation])
AC_ARG_WITH([root], AC_HELP_STRING([--with-root], [path to ROOT installation]), 
            [ROOTPATH=$with_root], [with_root=yes])
if test "x$ROOTPATH" = "xno"; then
    # When --without-root is given, do not use ROOT
    AC_MSG_NOTICE([Not using ROOT])
else
    ## Dumb default ROOT include and lib paths
    DEFAULTROOTINCPATH="/usr/include"
    DEFAULTROOTLIBPATH="/usr/lib"
    if test -d "$ROOTPATH"; then
        DEFAULTROOTINCPATH="$ROOTPATH/include"
        DEFAULTROOTLIBPATH="$ROOTPATH/lib"
    fi

    ## Use root-config to get better default paths
    AC_PATH_PROG(ROOTCONFIG, [root-config], [], [$ROOTPATH/bin:$PATH])
    if test -x "$ROOTCONFIG"; then
        DEFAULTROOTINCPATH=`$ROOTCONFIG --incdir`
        DEFAULTROOTLIBPATH=`$ROOTCONFIG --libdir`
    fi

    ## ROOT headers args
    AC_ARG_VAR(ROOTINCPATH,[path to directory containing ROOT headers @<:@ROOTPATH/include@:>@])
    AC_ARG_WITH([root-incpath],
      AC_HELP_STRING([--with-root-incpath],[path to directory containing ROOT headers @<:@ROOTPATH/include@:>@]),
      [ROOTINCPATH=$with_root_incpath], [ROOTINCPATH=${DEFAULTROOTINCPATH}])

    ## ROOT libraries args
    AC_ARG_VAR(ROOTLIBPATH,[path to directory containing ROOT libraries @<:@ROOTPATH/lib@:>@])
    AC_ARG_WITH([root-libpath],
      AC_HELP_STRING([--with-root-libpath],[path to directory containing ROOT libraries @<:@ROOTPATH/lib@:>@]),
      [ROOTLIBPATH=$with_root_libpath], [ROOTLIBPATH=${DEFAULTROOTLIBPATH}])

    ## Default ROOT build flags
    ROOTCPPFLAGS="-I$ROOTINCPATH"
    ROOTCXXFLAGS="-Wno-long-long"
    ROOTLDFLAGS=""
    ROOTLIBS="-L$ROOTLIBPATH -lCore -lCint -lRIO -lNet -lHist -lGraf -lGraf3d -lGpad -lTree"
    ROOTLIBS="$ROOTLIBS -lRint -lPostscript -lMatrix -lPhysics -lMathCore -lThread -pthread -lm -ldl -rdynamic"
    ## Use root-config to get specific flags for this version
    if test -x "$ROOTCONFIG"; then
        ROOTLIBS=`$ROOTCONFIG --libs`
        ROOTLDFLAGS=`$ROOTCONFIG --ldflags`
    fi

    ## Test ROOT headers
    AC_MSG_CHECKING([ROOT headers])
    if test -f ${ROOTINCPATH}/TROOT.h; then
      AC_MSG_RESULT([using ROOT headers in $ROOTINCPATH])
    else
      AC_MSG_WARN([ROOT headers not found under $ROOTINCPATH... disabling ROOT])
      with_root=no
    fi

    ## Test ROOT libs
    AC_MSG_CHECKING([ROOT libraries])
    if test -f ${ROOTLIBPATH}/libCore.so; then
      DEFAULTCXXFLAGS="-O2"
      AC_MSG_RESULT([using ROOT libraries in $ROOTLIBPATH])
    else
      AC_MSG_WARN([ROOT libraries not found under $ROOTLIBPATH... disabling ROOT])
      with_root=no
    fi

    AC_SUBST([ROOTPATH])
    AC_SUBST([ROOTCPPFLAGS])
    AC_SUBST([ROOTCXXFLAGS])
    AC_SUBST([ROOTLDFLAGS])
    AC_SUBST([ROOTLIBS])
fi

AM_CONDITIONAL(WITH_ROOT, [test x$with_root != xno])


################################


## Set default compiler flags if CXXFLAGS not defined
if test -z "$CXXFLAGS"; then
  CXXFLAGS=$DEFAULTCXXFLAGS
fi


AM_CXXFLAGS="-ansi -pedantic"
AC_SUBST(AM_CPPFLAGS)
AC_SUBST(AM_CXXFLAGS)
AC_SUBST(CXXFLAGS)

AC_CONFIG_FILES(Makefile)
AC_CONFIG_FILES(src/Makefile)
AC_CONFIG_FILES(include/Makefile)
AC_CONFIG_FILES(include/StatPatternRecognition/Makefile)

AC_OUTPUT
