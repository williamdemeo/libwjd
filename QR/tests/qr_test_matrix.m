function [A,Q,R,Rpiv] = qr_test_matrix(m, n, cnd, matrixfile, resultsfile)
% Create a random m by n matrix with singular values ranging from 
% 1 up to cnd (so that condition number of the matrix is cnd).
% Decompose the matrix using octave's qr function, both with and 
% without pivoting. Display the execution time.  Write the matrix 
% to a datafile called matrixfile.  Write the results of octave's 
% qr decomposition of the matrix to the file resultsfile.
%
% Inputs: 
%
%   m, n = numbers of rows, columns in test matrices
%          m should be at least n
%
%   cnd = (default=100) condition number of test matrices to generate
%         (ratio of largest to smallest singular value)
%         cnd should be at least 1
%
%   matrixfile   a string name of the file to which matrix will be written.
%   resultsfile  (optional) a string name of the file to which qr 
%                decomposition results will be written.  If no name is
%                given, results will not be written to a file.
%
% Outputs:
%   A  the matrix
%   Q  the Q in the QR decomposition A = QR.
%   R  the R in the QR decomposition A = QR.
%   Rpiv, E = the R and permutation matrix E from QR decomposition 
%             of A with pivoting
%
% Example: [A,Q,R,Rpiv] = qr_test_matrix(4,3,100,'matrix_4x3.txt', 'results_4x3.txt');
%
% Created 1997.11.28 by William J. DeMeo.
% Updated 2011.08.18.
%
write_results=true;
if nargin < 5,
  %resultsfile='resultsfile.txt';
  write_results=false;
  if nargin < 4,
    matrixfile='matrixfile.txt';
    if nargin < 3,
      cnd=100;
      if nargin < 2,
        n=8; m=8;
      end;
    end;
  end;
end;

if cnd<1, Error('Usage: test_matrix(m,n,cnd) with cnd not less than 1.'); end;

if m<n, Error('Usage: test_matrix(m,n,cnd) with m not less than n.'); end;

format long;

% Generate random matrix A, starting with the SVD of a random matrix
A=randn(m,n);    % Construct a random mxn matrix.
[u,s,v]=svd(A);  % Find its SVD.

% Transform the matrix so that the singular values range 
% from 1 to cnd, with uniformly distributed logarithms.
sd = [1, cnd, exp(rand(1,n-2)*log(cnd))];
s = diag(sd);
A=u(:,1:n)*s*v';

% Write the matrix A to datafile:
fid = fopen(matrixfile,'w');  % write file (discard previous contents)
fprintf(fid,'%f\n',A);
fclose(fid);

tt = cputime;
% Perform QR without pivoting:
[Q,R] = qr(A);
et = cputime - tt;

disp(sprintf('time for octave qr decomp (no-pivot): %7.4f sec',et));

tt = cputime;
% Perform QR with pivoting:
[Qpiv, Rpiv, P] = qr(A);
et = cputime - tt;

if write_results,
  % Now write original matrix A, and resulting decompositions
  %  Q, R, Qpiv, Rpiv, etc. to files in human readable form
  % for comparison with results of qr_test.c program:

  fid = fopen(resultsfile,'w');  % write file (discard previous contents)
  fprintf(fid,'\nA = \n'); 
  fclose(fid);
  write_matrix(A, resultsfile, 'a');
  write_matrix(A, 'A.txt', 'w');

  fid = fopen(resultsfile,'a');  % append file
  fprintf(fid,'\nQ = \n'); 
  fclose(fid);
  write_matrix(Q, resultsfile, 'a');
  write_matrix(Q, 'Q.txt', 'w');

  fid = fopen(resultsfile,'a');
  fprintf(fid,'\nR = \n'); 
  fclose(fid);
  write_matrix(R, resultsfile, 'a');
  write_matrix(R, 'R.txt', 'w');

  fid = fopen(resultsfile,'a');
  fprintf(fid,'\nQpiv = \n'); 
  fclose(fid);
  write_matrix(Qpiv, resultsfile, 'a');
  write_matrix(Qpiv, 'Qpiv.txt', 'w');

  fid = fopen(resultsfile,'a');
  fprintf(fid,'\nRpiv = \n'); 
  fclose(fid);
  write_matrix(Rpiv, resultsfile, 'a');
  write_matrix(Rpiv, 'Rpiv.txt', 'w');
end;

disp(sprintf('time for octave qr decomp (pivot): %7.4f sec',et));
%
% end test_matrix.m


