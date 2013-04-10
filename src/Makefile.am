LIBS = @LIBS@
LIBDIR = @LIBDIR@
LIBNAME = @LIBNAME@
my_include_paths = -I$(top_srcdir)/src -I$(top_srcdir) -I$(top_srcdir)/src/lr

AM_CXXFLAGS = $(CXXFLAGS) -I$(my_include_paths) $(BOOST_CPPFLAGS)
AM_CXXFLAGS += -DHAVE__BOOL  -D__NOPLUGIN__
AM_FFLAGS=$(my_include_paths) -x f77-cpp-input

my_ld_paths = $(LDFLAGS) -L$(top_srcdir)/src/$(LIBDIR) -L$(top_srcdir)/src/lr
my_ld_paths += -L/usr/local/lib
my_ld_flags = -lrjungle -llr -lgsl -lgslcblas

### SHARED
# the program to build
# bin_PROGRAMS = rjungle rjungle_static rjunglesparse rjunglesparse_static
# changed by bcw - 1/5/12 to only build shared binary
# building static binary again, but not sparse implementations
bin_PROGRAMS = rjungle rjungle_static

# compiler flags
# @rjungle_CPPFLAGS = $(my_include_paths)

# prog's sources
rjungle_SOURCES = main.cpp

# the local library needed to build the program (currently not installed
# in the lib directory, therefor relative path where it will be built!)
rjungle_LDADD = $(my_ld_paths) $(my_ld_flags) 

### STATIC
# compiler flags
# rjungle_static_CPPFLAGS = $(mycppflags)

# prog's sources
rjungle_static_SOURCES = main.cpp

# the local library needed to build the program (currently not installed
# in the lib directory, therefore relative path where it will be built!)
rjungle_static_LDADD = $(my_ld_paths) $(my_ld_flags) 

rjungle_static_LDFLAGS = -all-static

# additional include paths to support the library headers
INCLUDES = -I$(top_srcdir)/src -I$(top_srcdir)/src/$(LIBDIR) -I$(top_srcdir)/src/lr

# the subdirectories of src to go into
SUBDIRS = lr library
