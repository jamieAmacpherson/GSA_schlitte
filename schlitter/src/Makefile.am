#===============================================================================
# Makefile.am : for automake
#===============================================================================

bin_PROGRAMS = schlitter 

AM_CXXFLAGS = $(INTI_CFLAGS)

schlitter_SOURCES = \
	arg.c arg.h getgromacs.c getgromacs.h \
	safe.c safe.h schlitter.c schlitter.h

schlitter_LDADD = $(INTI_LIBS)

EXTRA_DIST = doxygen.cfg

CLEANFILES = $(TARGET) *.o 
DISTCLEANFILES = libtool config.cache config.log

