# Makefile for uniform_test

HOME = /accounts/grad/chip
UTIL = $(HOME)/lib/util/util_lib.a
INCL = -I$(HOME)/lib/include
MOMENT = $(HOME)/lib/moment
CC = cc

uniform_test: uniform_test.o unif.o $(MOMENT)/moment.o $(UTIL)
	$(CC) -o uniform_test uniform_test.o unif.o $(MOMENT)/moment.o $(UTIL) -lm

uniform_test.o: uniform_test.o
	$(CC) $(INCL) -c uniform_test.c 

unif.o: unif.c
	$(CC) $(INCL) -c unif.c 

$(MOMENT)/moment.o: $(MOMENT)/moment.c
	$(CC) -c -o $(MOMENT)/moment.o $(MOMENT)/moment.c

clean:
	rm *.o *~
