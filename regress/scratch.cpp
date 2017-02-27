/** reg() The main routine for running the regression and estimating B
 *                  in the equation XB = Y. N.B. this should be called only
 *                  \emph{after} you have already called qrpivot().
 *
 * @param M      number of rows of the design matrix X.
 * @param N      number of columns of X (n.b. we expect N < M).
 * @param QR     the QR decomposition of X, obtained by calling qrpivot().
 * @param leadu  on entry, the vector of leading u's resulting from qrpivot()
 *               on exit, the vector of coefficient estimates B, where y = XB.
 * @param E      the permutation matrix resulting from qrpivot().
 * @param y      an array of length M of "observables" (the rhs in XB = y).
 * @param B      on entry, an arbitrary length N vector
 *               on exit, the coefficient estimates.
 * @param cov    on entry, an arbitrary (pre-allocated) 1-d array of size NxN
 *               on exit, the covariance matrix.
 * @param se     on entry, an arbitrary (pre-allocated) array of length N
 *               on exit, the s.e.'s of the coefficient estimates.
 * @param e      on entry, an arbitrary (pre-allocated) array of length M
 *               on exit, the vector of residuals:  e = y - XB.
 * @param sigma  on exit, the mse =  y^te / (M-N).
 */
void Regression::reg(long M, long N, double *QR, double *leadu, double *E,
         double *y, double *B, double *cov, double *se, double *e, double *sigma){
}
