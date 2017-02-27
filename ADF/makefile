# Makefile for random

UTIL = ../util_lib.a
#NORMAL = $(HOME)/lib/random/normal
#UNIFORM = $(HOME)/lib/random/uniform
INCL = -I../include
#MOMENT = $(HOME)/lib/moment
CC = g++
COMPILE = -c 

dftest_test: ${UTIL} dftest_test.o dftest.o normal.o uniform.o moment.o regress.o ../QR/qr.o
	${CC} -o dftest_test  ${UTIL} dftest_test.o dftest.o normal.o uniform.o moment.o regress.o ../QR/qr.o -lblas -lm

%.o : %.c
	${CC} ${INCL} ${COMPILE} ${CFLAGS} $< -o $@

.PHONY: clean distclean

clean:
	rm -f *.o *.o~

distclean:
	rm -f *.o *.o~ ${APPS}



# random.o: random.c
# 	$(CC) $(INCL) -c random.c 

# normal.o: normal.c
# 	$(CC) $(INCL) -c -o normal.o normal.c 

# uniform.o: uniform.c
# 	$(CC) $(INCL) -c -o uniform.o uniform.c 

# moment.o: moment.c
# 	$(CC) -c -o moment.o moment.c

