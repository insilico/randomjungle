LIBS=@LIBS@
LIBDIR=@LIBDIR@
EXECNAME=@EXECNAME@
LIBNAME=@LIBNAME@
mycppflags=$(OPENMP_CXXFLAGS) $(BOOST_CPPFLAGS) -D__NOPLUGIN__

### SHARED
# the program to build
# bin_PROGRAMS = rjungle rjungle_static rjunglesparse rjunglesparse_static
# changed by bcw - 1/5/12 to only build shared binary
# building static binary again, but not sparse implementations
bin_PROGRAMS = rjungle rjungle_static

# compiler flags
@rjungle_CPPFLAGS = $(mycppflags)
#rjunglesparse_CPPFLAGS = $(mycppflags) -D__SPARSE_DATA__

# prog's sources
rjungle_SOURCES = main.cpp
#rjunglesparse_SOURCES = main.cpp

# the local library needed to build the program (currently not installed
# in the lib directory, therefor relative path where it will be built!)
rjungle_LDADD = $(OPENMP_LDFLAGS) -L$(top_srcdir)/src/$(LIBDIR) -L$(top_srcdir)/src/lr -lrjungle -llr @FLIBS@
#rjunglesparse_LDADD = $(OPENMP_CFLAGS) -L$(top_srcdir)/src/$(LIBDIR) -L$(top_srcdir)/src/lr -lrjunglesparse -llr @FLIBS@

### STATIC
# compiler flags
rjungle_static_CPPFLAGS = $(mycppflags)
#rjunglesparse_static_CPPFLAGS = $(mycppflags) -D__SPARSE_DATA__

# prog's sources
rjungle_static_SOURCES = main.cpp
#rjunglesparse_static_SOURCES = main.cpp

# the local library needed to build the program (currently not installed
# in the lib directory, therefor relative path where it will be built!)
rjungle_static_LDADD = -L$(top_srcdir)/src/$(LIBDIR) -L$(top_srcdir)/src/lr -lrjungle -llr @FLIBS@
#rjunglesparse_static_LDADD = -L$(top_srcdir)/src/$(LIBDIR) -L$(top_srcdir)/src/lr -lrjunglesparse -llr @FLIBS@

rjungle_static_LDFLAGS = -all-static $(OPENMP_CFLAGS)
#rjunglesparse_static_LDFLAGS = -all-static $(OPENMP_CFLAGS)

# additional include paths to support the library headers
INCLUDES = -I$(top_srcdir)/src -I$(top_srcdir)/src/$(LIBDIR) -I$(top_srcdir)/src/lr

# the subdirectories of src to go into
SUBDIRS = lr library
