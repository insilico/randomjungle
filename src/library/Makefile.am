LIBS=@LIBS@
LIBDIR=@LIBDIR@
LIBNAME=@LIBNAME@
myincludes = -I$(top_srcdir)/src -I$(top_srcdir) -I$(top_srcdir)/src/lr
mycppflags=$(OPENMP_CFLAGS) $(BOOST_CPPFLAGS) $(myincludes)
AM_FFLAGS=$(myincludes) -x f77-cpp-input

myheaders= \
BitMatrix.h BuildinGenerators.h ClassAtom.h CmpldTree.h CmplFct.h \
CondInf.h ErrorCodes.h IAM2WayImportance.h IAMClassAtom.h Importance.h JTreeTweaks.h \
TImportance.h RJungleTuneMtry.h RJungleOutlier.h RJungleProto.h \
DataFrame.h DataTreeSet.h Exception.h \
FittingFct.h Generator.h RJunglePrediction.h init.h\
SaveCollector.h PermImportance.h Prediction.h ioLines.h\
gzstream.h Helper.h INode.h JTreeCtrl.h librjungle.h Logistic.h Node.h \
Profiler.h Proximities.h RJungleAcc.h RJungleBinder.h RJungleCompiler.h \
RJungleConfusion.h RJungleCtrl.h RJungleFromXML.h RJungleGen.h \
RJungleGrow.h RJungleHelper.h RJungleImportance.h RJungleImpute.h \
RJungleIO.h RJunglePar.h RJungleProxi.h SClassAtom.h T2ClassAtom.h \
TermClassAtom.h TermResult.h TestClass.h TimeProf.h TMClassAtom.h \
LotusTermClassAtom.h treedefs.h Tree.h mvt.h

mysources= \
BitMatrix.cpp BuildinGenerators.cpp ClassAtom.cpp CmpldTree.cpp CmplFct.cpp \
CondInf.cpp IAM2WayImportance.cpp IAMClassAtom.cpp Importance.cpp JTreeTweaks.cpp \
TImportance.cpp RJungleTuneMtry.cpp RJungleOutlier.cpp RJungleProto.cpp \
DataFrame.cpp DataTreeSet.cpp Exception.cpp \
FittingFct.cpp RJunglePrediction.cpp\
SaveCollector.cpp PermImportance.cpp Prediction.cpp ioLines.cpp\
Generator.cpp gzstream.cpp Helper.cpp INode.cpp librjungle.cpp \
Logistic.cpp Node.cpp Profiler.cpp Proximities.cpp RJungleAcc.cpp RJungleBinder.cpp \
RJungleCompiler.cpp RJungleConfusion.cpp RJungleCtrl.cpp RJungleFromXML.cpp \
RJungleGen.cpp RJungleGrow.cpp RJungleHelper.cpp RJungleImportance.cpp \
RJungleImpute.cpp RJungleIO.cpp RJunglePar.cpp RJungleProxi.cpp \
SClassAtom.cpp T2ClassAtom.cpp TermClassAtom.cpp TermResult.cpp TestClass.cpp \
LotusTermClassAtom.cpp TimeProf.cpp TMClassAtom.cpp Tree.cpp mvt.cpp mvtfortran.f


# normal library build
#lib_LIBRARIES=librjungle.a librjunglesparse.a
#lib_LIBRARIES=librjungle.a

# compiler flags
#librjungle_a_CPPFLAGS=$(mycppflags)
#librjunglesparse_a_CPPFLAGS=$(mycppflags) -D__SPARSE_DATA__

# the list of header files that belong to the library
#librjungle_a_HEADERS=$(myheaders)
#ibrjunglesparse_a_HEADERS=$(myheaders)

# where to install the headers on the system
#librjungle_adir=$(includedir)/rjungle
#librjunglesparse_adir=$(includedir)/rjunglesparse

# the sources to add to the library and to add to the distribution
#librjungle_a_SOURCES=$(librjungle_a_HEADERS) $(mysources)
#librjunglesparse_a_SOURCES=$(librjunglesparse_a_HEADERS) $(mysources)

# shared library build with libtool
lib_LTLIBRARIES=librjungle-1.2.362.la
librjungle_1_2_362_la_SOURCES = $(mysources)

# compiler flags
librjungle_1_2_362_la_CPPFLAGS=$(mycppflags)

# the list of header files that belong to the library
librjungle_1_2_362_la_HEADERS=$(myheaders)

# where to install the headers on the system
librjungle_1_2_362_ladir=$(includedir)/rjungle
