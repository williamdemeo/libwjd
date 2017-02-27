# Master makefile for programs in the subdirectories under libwjd.
# Created on 20110819 by williamdemeo@gmail.com.
# Last modified on 20110819.
#

INCL = -Iinclude
LIBD = /usr/lib/atlas-base

CC = gcc
CFLAGS = -c -g
LDFLAGS = -g

qr_test: qr.o qr_test.o
	${CC} ${LDFLAGS} -o qr_test qr.o qr_test.o util_lib.a -L${LIBD} -lf77blas -lgfortran -latlas -lm

qr_test.o: qr_test.c
	${CC} ${CFLAGS} ${INCL} qr_test.c

qr.o: qr.c
	${CC} ${CFLAGS}  ${INCL} qr.c


qr.o: qr.c
	${CC} ${CFLAGS} ${INCL} -c qr.c

clean:
	rm -f *.o qr_test

#HDR = prototypes.h
# Additional dependencies
#qr.o : $(addprefix include/, ${HDR})
#qr_test.o : $(addprefix include/, ${HDR})

