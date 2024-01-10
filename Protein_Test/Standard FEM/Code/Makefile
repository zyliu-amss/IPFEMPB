PBE: PBE.o

PBE.o: PBE.c PBE.h Init.h func.h PBE_coef.h	

default: lib  nlnpoisson NPBE

all: lib   nlnpoisson

clean:
	-@/bin/rm -f PBE PBE.o

#include /opt/phg/phg-0.9.3/share/phg/Makefile.inc 
#include /opt/phg/phg-0.9.2/share/phg/Makefile.inc
#${PHG_MAKEFILE_INC}=/soft/apps/phg/intel-2017u4/impi-2017u3/phg-0.9.6/share/phg/Makefile.inc

#include ${PHG_MAKEFILE_INC}
include /soft/apps/phg/intel-2017u4/impi-2017u3/phg-0.9.6/share/phg/Makefile.inc
#include /soft/apps/phg/oneapi-2021.1.1/oneapimpi/phg-0.9.7-20230625-long_double-long/share/phg/Makefile.inc

.PHONY: default all clean distclean lib
