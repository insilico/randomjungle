### configure.ac                                        -*- Autoconf -*-
###
### Process this file with autoconf to produce a configure script.
###
### random jungle (RJ)
### Daniel Schwarz (c) 2008-2010
### http://www.randomjungle.com
###
### Content merged from different sources (autoscan, R, ...)
### RJ is free software; you can redistribute it and/or modify it under
### the terms of the GNU General Public License as published by the Free
### Software Foundation; either version 2 of the License, or (at your
### option) any later version.
###
### RJ is distributed in the hope that it will be useful, but WITHOUT ANY
### WARRANTY; without even the implied warranty of MERCHANTABILITY or
### FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public
### License for more details.
###

AC_PREREQ([2.60])

AC_INIT(randomjungle, [1.2.362])

AC_CONFIG_MACRO_DIR([m4])
AC_CONFIG_AUX_DIR([build-aux])
AC_CONFIG_SRCDIR(src/library/librjungle.cpp)
AC_CONFIG_HEADERS([config.h])

AM_INIT_AUTOMAKE

AC_SUBST(SRCDIR, src)
AC_SUBST(LIBDIR, library)
AC_SUBST(LIBLRDIR, lr)
AC_SUBST(MPISEP)
AC_SUBST(MPIINFIX)

AC_CANONICAL_BUILD
AC_CANONICAL_HOST

# Checks for programs and sets appropriate environment variables.
AC_PROG_CXX
AC_PROG_CPP
AC_PROG_LN_S
AC_PROG_LIBTOOL
AC_F77_LIBRARY_LDFLAGS

LT_INIT

# add parameter(s)
AC_ARG_ENABLE(debug, AC_HELP_STRING([--enable-debug],[Compile code for debugging.]), [debug_rjungle=yes])

# Set debug variables
if test "$debug_rjungle" == 'yes' ; then
    CPPFLAGS="-O0 -gdwarf-2 -Wall $CPPFLAGS -D__DEBUG__"
    CXXFLAGS="-O0 -gdwarf-2 -Wall $CXXFLAGS -D__DEBUG__"
fi

# Checks for header files.
AC_CHECK_HEADERS([math.h stdlib.h string.h sys/time.h getopt.h])

# Checks for typedefs, structures, and compiler characteristics.
AC_HEADER_STDBOOL
AC_C_INLINE
AC_TYPE_SIZE_T
AC_FUNC_ERROR_AT_LINE
AC_FUNC_STRTOD
AC_CHECK_FUNCS([round pow sqrt floor memset], [], AC_MSG_ERROR([required function missing]))
AC_CHECK_SIZEOF(void *)

# makeinfo
AC_CHECK_PROG([MAKEINFO], [makeinfo], [makeinfo], AC_MSG_ERROR([could not find required program: makeinfo]))

AC_MSG_NOTICE([--------------------------------------------------])

AC_CHECK_LIB(pthread, pthread_create)
LDFLAGS="$LDFLAGS -lpthread"

# Checks for libraries.
# zlib.m4
CHECK_ZLIB

# boost.m4
BOOST_REQUIRE
BOOST_FIND_HEADER(boost/dynamic_bitset.hpp, AC_MSG_ERROR([could not find required headers of Boost c++. Install boost c++ libraries.],))

# gsl.m4
AX_PATH_GSL(1.8,, AC_MSG_ERROR([could not find required version of GSL]))

# libxml2.m4
AM_PATH_XML2(2.0.0,, AC_MSG_ERROR([could not find required version of XML2]))

# Check for XML delegate library.
AC_ARG_WITH([xml],
    [AC_HELP_STRING([--without-xml],
                    [disable XML support])],
    [with_xml=$withval],
    [with_xml=$have_x])

if test "$with_xml" != 'yes' ; then
    DISTCHECK_CONFIG_FLAGS="${DISTCHECK_CONFIG_FLAGS} --with-xml=$with_xml "
fi

have_xml='no'
XML_LIBS=''
if test "$with_xml" != 'no'; then
    AC_MSG_CHECKING([for XML])
    AC_MSG_RESULT([])
    PERSIST_LDFLAGS=$LDFLAGS
    PERSIST_CPPFLAGS=$CPPFLAGS
    xml2_config=''
    AC_CHECK_PROGS(xml2_config,xml2-config,)dnl
    if test -n "$xml2_config"; then
        # Debian installs libxml headers under /usr/include/libxml2/libxml with
        # the shared library installed under /usr/lib, whereas the package
        # installs itself under $prefix/libxml and $prefix/lib.
        xml2_prefix=`xml2-config --prefix`
        if test -d "${xml2_prefix}/include/libxml2"; then
            CPPFLAGS="$CPPFLAGS -I${xml2_prefix}/include/libxml2"
        fi
        if test "${xml2_prefix}" != '/usr'; then
            LDFLAGS="$LDFLAGS -L${xml2_prefix}/lib"
        fi
    fi
    failed=0
    passed=0
    AC_CHECK_HEADER(libxml/parser.h,passed=`expr $passed + 1`,failed=`expr $failed + 1`)
    AC_CHECK_LIB(xml2,xmlParseExternalEntity,passed=`expr $passed + 1`,failed=`expr $failed + 1`,)
    AC_MSG_CHECKING([if XML package is complete])
    if test $passed -gt 0; then
        if test $failed -gt 0; then
            AC_MSG_RESULT([no -- some components failed test])
            have_xml='no (failed tests)'
            LDFLAGS="$PERSIST_LDFLAGS"
            CPPFLAGS="$PERSIST_CPPFLAGS"
        else
            XML_LIBS='-lxml2'
            LIBS="$XML_LIBS $LIBS"
            AC_DEFINE(XML_DELEGATE,1,Define if you have XML library)
            AC_MSG_RESULT([yes])
            have_xml='yes'
        fi
    else
        AC_MSG_RESULT([no])
    fi
fi
AM_CONDITIONAL(XML_DELEGATE,test "$have_xml" = 'yes')
AC_SUBST(XML_LIBS)

dnl ===========================================================================

# check for getopt in a separate library
AC_DEFINE([HAVE_GETOPT_LONG], [], [Defined when getopt_long is available])
AC_CHECK_FUNC([getopt_long], AC_DEFINE([HAVE_GETOPT_LONG]),)

#HAVE_GETOPT=no
#AC_CHECK_LIB(getopt, getopt, HAVE_GETOPT="yes")
#if test "x$HAVE_GETOPT" = xyes ; then
#  GETOPT_LIBS="-lgetopt"
#  AC_SUBST(GETOPT_LIBS)
#fi
#
#dnl Check for getopt_long (Solaris)
#AC_CHECK_FUNCS([getopt_long],,AC_CHECK_LIB([iberty],[getopt_long],OTHERLIBS="$OTHERLIBS -liberty"))
#AC_SUBST(OTHERLIBS)
#
#AC_DEFINE([HAVE_GETOPT_LONG], [], [Defined when getopt_long is available])
#AC_CHECK_FUNC([getopt_long], AC_DEFINE([HAVE_GETOPT_LONG]),
#  AC_CHECK_LIB([gnugetopt],[getopt_long], ,
#    AC_MSG_ERROR([*** getopt_long not available!])))


AC_MSG_NOTICE([--------------------------------------------------])

# Process makefiles
AC_CONFIG_FILES([Makefile doc/Makefile src/lr/Makefile src/library/Makefile src/Makefile])
AC_OUTPUT
