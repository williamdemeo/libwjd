# makefile for QR and QRpiv
# makefile for ctimedgemm.cc
# modified: 20110816, 20031123
#
# HOME = /home/william/pub/DeMeo/computing/src/C
HOME = /home/williamdemeo/git/PROGRAMMING/C
UTIL = $(HOME)/libwjd/util_lib.a
# INCL = -I$(HOME)/libwjd/include -I$(HOME)/stat243/samples
INCL = -I$(HOME)/libwjd/include
# LIBD = /usr/local/atlas/lib
LIBD = /usr/lib/atlas-base
# CC = cc
DECOMP = $(HOME)/libwjd/QR
CC = gcc

# Variables for optimization on SUN
# Reference:  Sunsoft Performace Library
# /opt/SUNWspro/READMEs/libsunperf_preface.ps (p. viii)
#LINK = -lsunperf -dalign -lblas # use these for link command
COMPILE = -dalign # use this option for compile command on single processor
MPCOMPILE = -dalign -mt  # multiprocessor compile command

# 2001.10.21: No longer using the Sunsoft Performance Library
#             Define new LINK variable here:
# Defines for use with ATLAS are in Make.inc
# include Make.inc
QRpiv: QRpiv.o $(DECOMP)/House.o
	$(CC) -o QRpiv QRpiv.o $(DECOMP)/House.o $(UTIL) -L$(LIBdir) -lf77blas -lg2c -latlas -lm

QRpiv.o: QRpiv.c
	$(CC) $(COMPILE) $(INCL) -c QRpiv.c

QR: QR.o $(DECOMP)/House.o
	$(CC) -o QR QR.o $(DECOMP)/House.o $(UTIL) -L$(LIBdir) -lf77blas -lg2c -latlas -lm

QR.o: QR.c
	$(CC) $(COMPILE) $(INCL) -c QR.c

QRptime: QRptime.o $(DECOMP)/House.o
	$(CC) -o QRptime QRptime.o $(DECOMP)/House.o $(UTIL) -L$(LIBdir) -lf77blas -lg2c -latlas -lm

QRptime.o: QRptime.c
	$(CC) $(COMPILE) $(INCL) -c QRptime.c

QRtime: QRtime.o $(DECOMP)/House.o
	$(CC) -o QRtime QRtime.o $(DECOMP)/House.o $(UTIL) -L$(LIBdir) -lf77blas -lg2c -latlas -lm

QRtime.o: QRtime.c
	$(CC) $(COMPILE) $(INCL) -c QRtime.c

$(DECOMP)/House.o: $(DECOMP)/House.c
	$(CC) $(COMPILE) $(INCL) -c -o $(DECOMP)/House.o $(DECOMP)/House.c

clean:
	rm -f *.o *.a *~ 


