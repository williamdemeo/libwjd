# makefile for QR and QRpiv
# makefile for ctimedgemm.cc
# modified: 20110816, 20031123
#
# Based on GNU/Linux Application Programming 2ed, by Tim Jones, p. 68.

SRC_FILES=QR_test.c House.c

VPATH = ../QR
LIBD = /usr/lib/atlas-base
# LIBD = /usr/local/atlas/lib
CC = gcc
CFLAGS = -c -g
LDFLAGS = -g

qr_test: ${OBJ_FILES}
	${CC} ${LDFLAGS} -o qr_test ${OBJ_FILES} -L${LIBD} -lf77blas -lgfortran -latlas -lm

%.o:%.c
	${CC} ${CFLAGS} -o $@ $

clean:
	rm -f *.o qr_test
