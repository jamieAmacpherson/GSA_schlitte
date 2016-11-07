#===============================================================================
# Makefile.am : for automake
#===============================================================================

bin_PROGRAMS = schlitter 

AM_CXXFLAGS = $(INTI_CFLAGS)

schlitter_SOURCES = \
arg.c arg.h argpdb.h config.h error.c error.h getpdb.c getpdb.h \
gettraj.c gettraj.h modstring.c modstring.h pattern.c pattern.h \
pdb_structure.h safe.c safe.h schlitter.c schlitter.h \
seq.h vector.c vector.h

schlitter_LDADD = $(INTI_LIBS)

EXTRA_DIST = doxygen.cfg

CLEANFILES = $(TARGET) *.o 
DISTCLEANFILES = libtool config.cache config.log

