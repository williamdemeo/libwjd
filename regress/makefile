# Makefile for regress.c
# Created on 2011.09.22 by williamdemeo@gmail.com
UTIL = ../util_lib.a
DECOMP = ../QR
INCL = -I../include
#LIBD = /usr/local/ATLAS/lib
LIBD = 
#LIBS = -lf77blas -lcblas -lgfortran -latlas -lm
LIBS = -llapack -lblas
CC = gcc
CFLAGS = -c -g
LDFLAGS = -g
COMPILE = -O3

regress: regress.o ${DECOMP}/qr.o
	${CC} ${LDFLAGS} -o regress regress.o ${DECOMP}/qr.o ${UTIL} -L${LIBD} ${LIBS}
#	${CC} ${LDFLAGS} -o regress regress.o ${DECOMP}/qr.o ${UTIL} -L${LIBD} -L${LIBL} -L${LF2C} ${LIBS}

${DECOMP}/qr.o: ${DECOMP}/qr.c
	${CC} ${COMPILE} ${INCL} -c -o ${DECOMP}/qr.o ${DECOMP}/qr.c

%.o : %.c
	${CC} ${COMPILE} ${CFLAGS} ${INCL} $< -o $@

.PHONY: clean distclean

clean:
	rm -f *.o

distclean:
	rm -f *.o ${APPS}


# Additional dependencies
H_DIR := ../include
H_FILES := $(addprefix $(H_DIR)/,prototypes.h timer.h)
$(DECOMP)/qr.o : ${H_FILES}
regress.o : ${H_FILES}



#LIBL = /usr/local/CLAPACK-3.2.1
#LF2C = /usr/local/CLAPACK-3.2.1/F2CLIBS
#LIBS = -llapack -lcblas -lf77blas -lf2c -lgfortran  -latlas -lm
#LIBS = -llapack_LINUX.a -lcblaswr -lcblas -lf77blas -lf2c -latlas -lm
#LIBS = /usr/local/CLAPACK-3.2.1/lapack_LINUX.a -lcblas -lf2c -latlas -lf77blas -lgfortran -lm
#LIBS = -lcblaswr -lcblas -lf2c -lf77blas -lgfortran -latlas -lm
#LIBS = -lf77blas -lgfortran -latlas -lm
