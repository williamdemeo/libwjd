/*  bits.c  utility for manipulating bits  

Created by William DeMeo
on Frebruary 26, 1999
Last modified: February 26, 1999
*/
#include<limits.h>
#include<stdlib.h>

/* powers of two  */
unsigned int Powers[] = {1, 2, 4, 8, 16, 32, 64, 128};

void bit_print(int a);

/* for using bit_print on its own */
/*
void 
main(int argc, char *argv[]) {
     int instr;
     if (argc < 2 || argc > 3) {
	  printf("Wrong number of arguments\n");
	  exit(1);
     }

     instr = atoi(argv[1]); 
     printf("Decimal: %d\n",instr);
     printf("Binary:  ");
     bit_print(instr);
     printf("\n");
}
*/
void bit_print (int a) {
     int i;
     int n = sizeof(int) * CHAR_BIT;  /* defined in limits.h */
     int mask = 1 << (n-1);           /* mask = 1000...0 */

     for(i=1;i<=n;++i){
	  putchar(((a & mask) == 0) ? '0' : '1');
	  a <<= 1;  /* shift a one bit */
	  if (i % CHAR_BIT == 0 && i < n) putchar(' ');
     }
}
