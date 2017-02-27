# Makefile for QR programs.
# Created on 20110819 by williamdemeo@gmail.com.
# Last modified on 20110819.
#
# Remarks: as of 20110819 it's (finally) working properly.
#
APPS = qr_test qrpivot_test qr_time qrpivot_time
#INCL = -I../include -I/usr/local/atlas-3.8.4/include
INCL = -I../include
LIBD = 
#LIBD = /usr/local/ATLAS/lib
#LIBS = -lf77blas -lcblas -lgfortran -latlas -lm
LIBS = -llapack -lblas 

## for C use these:
# UTIL = ../util_lib.a 
# CC = gcc

## for C++ use these:
UTIL = ../util_lib_cc.a 
CC = g++

CFLAGS = -c -g
LDFLAGS = -g
COMPILE = -O3

all : ${APPS}

qr_test: qr.o qr_test.o ${UTIL}
	${CC} ${LDFLAGS} -o qr_test qr.o qr_test.o ${UTIL} -L${LIBD} ${LIBS}

qr_time: qr.o qr_time.o
	${CC} ${LDFLAGS} -o qr_time qr.o qr_time.o ${UTIL} -L${LIBD} ${LIBS}

qrpivot_test: qr.o qrpivot_test.o
	${CC} ${LDFLAGS} -o qrpivot_test qr.o qrpivot_test.o ${UTIL} -L${LIBD} ${LIBS}

qrpivot_time: qr.o qrpivot_time.o
	${CC} ${LDFLAGS} -o qrpivot_time qr.o qrpivot_time.o ${UTIL} -L${LIBD} ${LIBS}

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
qr.o : ${H_FILES}
qr_test.o : ${H_FILES}
qr_time.o : ${H_FILES}
qrpivot_test.o : ${H_FILES}



#LIBD = /usr/local/CLAPACK
#LIBF2C = /usr/local/CLAPACK/F2CLIBS

#LIBS = -lcblas -latlas /usr/local/CLAPACK/F2CLIBS/libf2c.a -lm
#LIBS = ${LIBD}/lapack_LINUX.a ${LIBD}/F2CLIBS/libf2c.a -lm
#LIBS = ${LIBD}/lapack_LINUX.a -lgfortran -lm
#LIBS = ${LIBD}/lapack_LINUX.a -lcblas ${LIBD}/F2CLIBS/libf2c.a -lm
#LF2C = /usr/local/CLAPACK-3.2.1/F2CLIBS
#LIBL = /usr/local/CLAPACK-3.2.1
#UTIL = ../util_lib_cc.a 
#LIBS = -llapack -lcblas -lf77blas -lgfortran -latlas -lm
#LIBS = -lcblaswr -lcblas -lf2c -lf77blas -lgfortran -latlas -lm
#LIBS = /usr/local/CLAPACK-3.2.1/lapack_LINUX.a -lcblas -lf2c -latlas -lf77blas -lgfortran -lm
#LIBS = -llapack -lcblas -lf77blas -lf2c -latlas -lm
#LIBS = -llapack_LINUX -lcblaswr -lcblas -lf77blas -lf2c -latlas -lm
#LIBS = -llapack -lcblas -lf77blas -lf2c -latlas -lm
#LIBS = -lf77blas -lgfortran -latlas -lm
