Notes on testing the QR subroutines.
williamdemeo@gmail.com
2011.08.19

The following assumes you're in the QR directory which contains
a subdirectory called tests.

First create the test matrices in octave, as follows:

Change to the tests directory:

$ cd ./tests

Start octave:

$ octave

octave:> qr_test_matrix(8,8);

This constructs an 8x8 matrix with condition number 
100 (which can be changed) and stores it in the file called
matrixfile.txt.  It then uses octave's QR routine to compute
the QR decomposition and stores the results in small_matrix_QR.txt

Now change back to the QR directory and run

$ ./qr_test

The program will ask for the filename (e.g. ./tests/matrixfile.txt) and the 
dimensions of the matrix (e.g. 8 x 8).  It outputs ot the terminal the 
results of the qr() routine in qr.c.

A brief review of condition numbers is given below in the 
section CONDITION NUMBERS.


For timing purposes, you could repeat the foregoing with much
larger matrices, say, 512x512, instead of 8x8.  E.g.,

octave:> qr_test_matrix(512,512,100,'bigmatrixfile.txt','bigresultsfile.txt');

$ ./qr_time

Currently, this routine doesn't prompt for input file or matrix dimensions.  
They are hard-coded as 'bigmatrixfile.txt and 512x512.

(On my laptop, timing tests showed:

octave takes 0.45 sec to decompose a 512x512 matrix (no pivot).
MATLAB takes 0.15 sec to decompose the same matrix (no pivot).
qr.c   takes 0.17 sec to decompose the same matrix (no pivot).



=================
CONDITION NUMBERS
-----------------
Recall, the condition number of Ax=b is the ratio
of the largest/smallest singular values of A.
That is, the ratio of the largest/smallest evalues of 
the square matrix A'A.  

To solve the system Ax=b, we form the normal equations
A'Ax=A'b, and then invert the matrix A'A.
If the condition number is large, this means the columns
of A are nearly linearly independent and the solution

        x = (A'A)^{-1}A'b

will not be robust in the sense that small changes in b 
may lead to large changes in x.

See also: [1] and [2]



[1] http://en.wikipedia.org/wiki/Condition_number
[2] Demmel, "Applied Numerical Linear Algebra," page 4.