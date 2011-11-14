LIBS=@LIBS@
LIBDIR=@LIBDIR@
LIBNAME=@LIBNAME@
myincludes = -I$(top_srcdir)/src -I$(top_srcdir)
mycppflags=$(OPENMP_CFLAGS) $(myincludes)

myheaders= \
amar.h ambs.h amdym.h amdyv_array.h amdyv.h amiv.h amma.h am_string_array.h am_string.h amvpv.h file.h genarray.h lin_conjgrad.h lr.h lrutils.h pfload.h prec.h precs.h score.h spardat.h standard.h

mysources= \
amar.c ambs.c amdym.c amdyv_array.c amdyv.c amiv.c am_string_array.c am_string.c amvpv.c file.c genarray.c lin_conjgrad.c lr.c lrutils.c pfload.c prec.c precs.c score.c spardat.c

# normal library build
lib_LIBRARIES=liblr.a

# compiler flags
#liblr_CPPFLAGS=$(mycppflags)

# the list of header files that belong to the library
liblr_a_HEADERS=$(myheaders)

# where to install the headers on the system
liblr_adir=$(includedir)/rjungle

# the sources to add to the library and to add to the distribution
liblr_a_SOURCES=$(liblr_a_HEADERS) $(mysources) 

