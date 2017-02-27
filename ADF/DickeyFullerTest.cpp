/**
 * @file
 * @author William
 *
 * Copyright 2011 Shooting Soul, LLC. and William DeMeo  All rights reserved.
 */

#include <iostream>
#include <string>
#include <sstream>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <limits>
#include "f2c.h"
#include "clapack.h"
#include "prototypes.h"
#include "tools.h"
#include "DickeyFullerTest.h"
#include "RandomSequence.h"
#include "UnitTest.h"

//#define DEBUG_DickeyFullerTest 1

#ifdef DEBUG_DickeyFullerTest
	#include <assert.h>
#endif

namespace ttm {

DickeyFullerTest::DickeyFullerTest(size_t N, double* y) :
		m_length(N), m_y(y), m_model(MODEL_AR), m_lags(1), m_alpha(0.05){


	long n=(long)m_length, INC = (long)1;
	double* diffy = new double[m_length];  // diffy = dmalloc(N);
	dcopy_(&n, m_y, &INC, diffy, &INC);  // diffy <- m_y

	// Get first-differenced series diffy(0) = y1-y0, diffy(1) = y2-y1, etc
	difference(m_length, diffy);
	m_diffy=diffy;
	// m_diffy now contains the first differences of m_y.  It has length  m_legnth-1.

	m_effective_sample_size = m_length - m_lags;  // effective sample size accounting for lagged differences

}


DickeyFullerTest::DickeyFullerTest(size_t N, double* y, model_type m, size_t l, double a) :
				m_length(N), m_y(y), m_model(m), m_lags(l), m_alpha(a) {

	long n=(long)m_length, INC = (long)1;
	double* diffy = new double[m_length];  // diffy = dmalloc(N);
	dcopy_(&n, m_y, &INC, diffy, &INC);  // diffy <- m_y

	// Get first-differenced series diffy(0) = y1-y0, diffy(1) = y2-y1, etc
	difference(m_length, diffy);
	m_diffy=diffy;
	// m_diffy now contains the first differences of m_y.  It has length  m_legnth-1.

	m_effective_sample_size = m_length - m_lags;   // effective sample size accounting for lagged differences
}

DickeyFullerTest::~DickeyFullerTest() {

}

const string DickeyFullerTest::MODEL_NAME[] = { "AR", "ARD" };

// The next three are defined already in .h file.
//void DickeyFullerTest::set_model(model_type m){	m_model = m; }
//void DickeyFullerTest::set_lags(size_t l){	m_lags = l; }
//int DickeyFullerTest::get_test_H(){	return m_test_H; }

//void DickeyFullerTest::display_regression_results() const{
//	display_regression_results(std::cout);
//}

void DickeyFullerTest::display_regression_results(ostream *out) const{
	ttm_assert(m_reg!=NULL,"m_reg not initialized");

	m_reg->display_results(out);
}

void DickeyFullerTest::check_initial_member_variables() const{
	if (m_length==NULL
			|| m_y==NULL
			|| m_diffy==NULL
		//	|| m_model==NULL
			|| m_lags==NULL
			|| m_effective_sample_size==NULL
			|| m_alpha==NULL)
	{
		ttm_throw ("some essential member variables not initialized")
	}
}

// These should be non-null after running test.
void DickeyFullerTest::check_all_member_variables() const{

	check_initial_member_variables();

	if (m_test_H==NULL || m_p_value==NULL || m_test_statistic==NULL || m_reg==NULL)
	{
		ttm_throw ("some member variables not initialized")
	}
}


//void DickeyFullerTest::display_attributes() const{
//	display_attributes(&std::cout);
//}

void DickeyFullerTest::display_attributes(ostream *out) const{
	check_initial_member_variables();
	*out << "\n DickeyFullerTest attributes:" << endl;
	*out << "\t m_length: \t\t "<< m_length << endl;
	*out << "\t m_model: \t\t" << MODEL_NAME[m_model] << endl;
	*out << "\t m_lags: \t\t" << m_lags << endl;
	*out << "\t m_effective_sample_size: \t" << m_effective_sample_size << endl;
	*out << "\t m_alpha: \t\t" << m_alpha << endl;
}

//void DickeyFullerTest::display_results() const{
//	display_results(&std::cout);
//}

void DickeyFullerTest::display_results(ostream *out) const{
	check_all_member_variables();
	*out << "\n DickeyFullerTest results:" << endl;
	*out << "\t m_test_H: \t\t" << m_test_H << endl;
	*out << "\t m_p_value: \t\t" << m_p_value << endl;
	*out << "\t m_test_statistic: \t\t" << m_test_statistic << endl;
	*out << "\t regression results:" << endl;
	m_reg->display_results(out);

}

void DickeyFullerTest::copy_first_differences(double* diffy) const{
	if(diffy==NULL) {
		ttm_throw("invalid input");
	}
	if(m_length==NULL || m_diffy==NULL){
		ttm_throw("m_length or m_diffy not initialized");
	}
	long n=(long)m_length-1, INC = (long)1;
	dcopy_(&n, m_diffy, &INC, diffy, &INC);  // diffy <- m_diffy

}

// Test the statistic
// run_test(const double& m_alpha, const double& testStat, const long& sampleSize, const  model_type& Model, int& testH, double& pVal)
void DickeyFullerTest::run_test(){
#ifdef DEBUG_DickeyFullerTest
	std::cout << "\t run_test()" << endl;
#endif
  int row, col;
  double testCValue=-999; // initialize to -999 so we don't accidentally reject the null if testCValue doesn't get asigned a value.

  compute_test_statistic();  // Recompute m_test_statistic each time we call run_test()
                             // to be sure we're not using an old m_test_statistic
                             // associated with different modeling assumptions.

  // upon return of this fn, row is the row of the CVTable we need:
  m_CV_sample_size = get_sample_size_index(m_effective_sample_size, row);

  //  upon return of this fn, col is the col of the CVTable we need:
  get_sig_index(m_alpha, col);

  switch(m_model){
  case MODEL_AR:
	  // CVTable's have size 15x50=CVTableNumRows*CVTableNumCols...
	  m_c_value = CVTable_AR[row*CVTableNumCols + col]; // ...stored in ROW-major format!
	  compute_p_value(CVTable_AR+(row*CVTableNumCols));
	  break;
  case MODEL_ARD:
	  m_c_value = CVTable_ARD[row*CVTableNumCols + col];
	  compute_p_value(CVTable_ARD+(row*CVTableNumCols));
	  break;
  default:
	  // TODO(wjd): fix this!!
	  //ttm_throw(string("unkonwn model type: ")+string((int)m_model));
	  ttm_throw("unkonwn model type");
	  break;
  }


  if (m_test_statistic < m_c_value)
    m_test_H = 1;  // reject the null in favor of the alternative
  else
    m_test_H = 0;  // fail to reject the null

}


// Typically we won't call this routine directly.  It will be called by run_test().
void DickeyFullerTest::compute_test_statistic(){
#ifdef DEBUG_DickeyFullerTest
	std::cout << "\t compute_test_statistic()" << endl;
#endif

	double *x;
	long i, N, nrow, ncol, mindim, INC=1;  // n=N since blas doesn't like const argument.
	int extra_cols;   // how many extra columns to add because of drift term and/or time-trend term
	                  // for simple AR model extra_cols=0
				      // for ARD model extra_cols=1
                      // Add 1 more to the above values of extra_cols, when there is a time-trend term.
                      // In each case, the coefficient on y_{t-1} in the regression will be m_B[extra_cols];

	check_initial_member_variables();

  if (m_model==MODEL_AR){  extra_cols=0; }
  else if(m_model==MODEL_ARD){ extra_cols=1; }
  else{
	  ttm_throw("unkonwn model type");
  }

  nrow = m_length-1;     // Effective sample size (after taking first differences)
  ncol = 1+extra_cols;        // if m_model==1, first column is a vector of 1's (for drift term)
  dassert(ncol<nrow,"ncol>=nrow not supported in this version");
  mindim = ncol;
  N = m_length;       // need a non-private value, to pass to blas routines and such.

  x = new double[nrow*(ncol+1)];  // x = dmalloc(nrow*(ncol+1));
  if (extra_cols==1){
	  // first column of x is a vector of 1's
	  for (i=0; i<nrow; ++i) x[i] = 1.0;  // storing elements in COLUMN-major format
  }
  dcopy_(&nrow, m_y, &INC, x+(nrow*extra_cols), &INC);    // x(:,1) <-  m_y(0:N-2)
  dcopy_(&nrow, m_diffy, &INC, x+(nrow*ncol), &INC);        // x(:,2) <- diffy(0:N-2) (right-hand-side)
  // x is the design matrix, with the independent variable (diffy) in the last column.
  /* That is, x is the matrix:      1    y0    y1-y0
                                    1    y1    y2-y1
                                    1    y2    y3-y2
				                      ...             ...in case m_model=1.
				                                      If m_model=0, there is no column of 1's.
				                                                                             */
#ifdef DEBUG_DickeyFullerTest
  std::cout << "Displaying design matrix [X, y] = " << endl;
  matprint(x,nrow, ncol+1);
#endif

  /// RUN REGRESSION ///
  //Regression m_reg((size_t)nrow, (size_t)ncol, x, m_diffy);
  m_reg = new Regression((size_t)nrow, (size_t)ncol, x, m_diffy);

#ifdef DEBUG_DickeyFullerTest
  std::cout << "Displaying design matrix X = " << endl;
  m_reg->display_design_matrix();
  //cout << "\nThe independent variable is:" << endl;	matprint(m_reg.m_y, m_nrow, 1);
  m_reg->display_results(&std::cout);
#endif

  m_reg->run_regression();   // might consider making it run automatically in constructor

  // Compute the test statistic:
  m_test_statistic = ( m_reg->get_coef(extra_cols) )/( m_reg->get_se(extra_cols) );
  // Later make this a function to accommodate the various kinds of tests we could do.

  delete[] x;
}



bool DickeyFullerTest::compute_p_value(const double *table){
#ifdef DEBUG_DickeyFullerTest
	std::cout << "\t compute_p_value()" << endl;
#endif

  int j;
  double fd;
  if (m_test_statistic<=table[0]){
    m_p_value = sigLevel[0];
    if (m_test_statistic<table[0]){
        std::cout << "\n\t WARNING: Test statistic, " << m_test_statistic << ", is outside of table values (p-value only approximated)." << endl;
    }
    return true;
  }
  j=0;
  while((m_test_statistic>table[j]) && (j<CVTableNumCols)){
    j++;
  }
  if (j==CVTableNumCols){
	    std::cout << "\n\t WARNING: Test statistic, " << m_test_statistic << ", is outside of table values (p-value only approximated)." << endl;
    m_p_value = sigLevel[j-1];
    return true;
  }
  // At this point m_test_statistic must be somewhere between the values table[j-1] and table[j].
  // Check that:
  if (m_test_statistic < table[j-1] || m_test_statistic > table[j]){
	  ttm_throw("unexpected test statistic\n p-value not set");
  }

  // Compute the "fractional distance" for the purpose of interpolating sig values:
  fd = (m_test_statistic - table[j-1])/(table[j] - table[j-1]);
  // When fd closer to 1 (0), m_test_statistic is closer to table[j] (table[j-1]).

  // Linear interpolation (to get approximate p-value).
  m_p_value = (1-fd)*sigLevel[j-1] + fd*sigLevel[j];
  // To get the most conservative p-value estimate, we could just take pVal = table[j].
  // In that case, if we reported a p-value of 0.020 = table[j], it would m_mean the
  // true p-value is somewhere in the interval (0.015, 0.020) = (table[j-1], table[j]).
  // This seems overly conservative, and interpolation is probably better.
  return(1);
}


// get_sig_index -- get index of the element in DFTest_sigLevel array (see dftest.h)
//                that is the greatest value not above m_alpha.
int DickeyFullerTest::get_sig_index(const double& alpha, int& i){
#ifdef DEBUG_DickeyFullerTest
	std::cout << "\t get_sig_index()" << endl;
#endif

  int j;
  i = -1;

#ifdef DEBUG_DickeyFullerTest
  assert(alpha>0.0009); assert(alpha<0.9999);
#endif

  if (alpha==sigLevel[0]){
    i=0;
    m_sig_level=alpha;
    return(i);
  }
  if (sigLevel[0]<alpha){
    j=0;
    while (sigLevel[j]<alpha){
      j++;
      dassert(j<sigLevelLength,"couldn't find appropriate significance level");
    }
    if (alpha < sigLevel[j]){
    	i=j-1;
       std::cout << "\n\t WARNING: Testing at significance level " << sigLevel[i] << " < alpha = " << alpha  << "... hope that's okay!" << endl;
    } else if (alpha == sigLevel[j]) {
    	i=j;
    } else {
    	ttm_throw("couldn't find appropriate significance level");
    }
    m_sig_level = sigLevel[i];
    return(i);
  }
  // If you made it to this point, something's wrong.
  ttm_throw("couldn't find appropriate significance level");
}



// get_sample_size_index -- get the index of the element in DFTest_CVsampleSize array (see dftest.h)
//                       that is the greatest value not above sampleSize.
int DickeyFullerTest::get_sample_size_index(const int& sampleSize, int &i){
#ifdef DEBUG_DickeyFullerTest
	std::cout << "\t get_sample_size_index()" << endl;
#endif

  int j;
  i = -1;
  dassert(sampleSize>9,"sample sizes less than 10 not supported");
  if (sampleSize==CVSampleSize[0]){   // sample size is 10
    i=0;
    return(CVSampleSize[0]);
  }
  if (CVSampleSize[0]<sampleSize){
    j=0;
    while (CVSampleSize[j]<sampleSize){
      j++;
      dassert(j<CVSampleSizeLength,"index out of bounds");
    }
    i=j; // i is the row of the CVTable we will use.

    // DFTest_CVSampleSize[i] should be the greatest sample size that is no greater than sampleSize. So,
    // i=j-1
    // On second thought, I think DFTest_CVSampleSize[i] should be the SMALLEST sample size that is
    // no SMALLER than than sampleSize. So,
    // i=j
    // This will result in more conservative test results.

    // Check that it's what we expect:
    if (sampleSize < CVSampleSize[i-1] || sampleSize > CVSampleSize[i]) {
      std::cout << "\n\t WARNING: couldn't find appropriate sample size for test...." << endl;
      std::cout << "\n\t WARNING:     ...sample size index may not be set appropriately." << endl;
    return(-1);
    } else {
      return(CVSampleSize[i]);
    }
  }
  // If you made it to this point, something's wrong.
  ttm_throw("couldn't find appropriate sample size");

}


const double DickeyFullerTest::sigLevel[] = { 0.001, 0.005, 0.010, 0.015, 0.020, 0.025, 0.030, 0.035, 0.040, 0.045, 0.050, 0.055, 0.060, 0.065, 0.070,
			     0.075, 0.080, 0.085, 0.090, 0.095, 0.100, 0.125, 0.150, 0.175, 0.200, 0.800, 0.825, 0.850, 0.875, 0.900,
			     0.905, 0.910, 0.915, 0.920, 0.925, 0.930, 0.935, 0.940, 0.945, 0.950, 0.955, 0.960, 0.965, 0.970, 0.975,
			     0.980, 0.985, 0.990, 0.995, 0.999};
const int DickeyFullerTest::sigLevelLength = 50;

const int DickeyFullerTest::CVSampleSize[] = {10, 15, 20, 25, 30, 40, 50, 75, 100, 150, 200, 300, 500, 1000, 10000};
const int DickeyFullerTest::CVSampleSizeLength= 15;
const int DickeyFullerTest::CVTableNumCols = sigLevelLength;
const int DickeyFullerTest::CVTableNumRows = CVSampleSizeLength;
const int DickeyFullerTest::CVTableLength = CVTableNumCols* CVTableNumRows;


// Note: the table is stored in ROW-major format!
const double DickeyFullerTest::CVTable_AR[] = {
  -3.9976, -3.1423, -2.7891, -2.5834, -2.4351, -2.3201, -2.2260, -2.1460, -2.0766, -2.0148, -1.9594, -1.9095, -1.8629, -1.8201, -1.7799, -1.7427, -1.7075,
  -1.6744, -1.6427, -1.6126, -1.5839, -1.4569, -1.3501, -1.2567, -1.1730,  0.4848,  0.5897,  0.7056,  0.8373,  0.9913,  1.0258,  1.0621,  1.0996,  1.1392,
   1.1810,  1.2254,  1.2722,  1.3230,  1.3775,  1.4361,  1.5004,  1.5720,  1.6520,  1.7437,  1.8514,  1.9818,  2.1485,  2.3836,  2.7840,  3.7494, //10

  -3.7444, -3.0262, -2.7138, -2.5289, -2.3954, -2.2909, -2.2044, -2.1305, -2.0659, -2.0083, -1.9560, -1.9084, -1.8648, -1.8239, -1.7858, -1.7501, -1.7161,
  -1.6841, -1.6537, -1.6248, -1.5970, -1.4738, -1.3690, -1.2773, -1.1950,  0.4548,  0.5589,  0.6733,  0.8024,  0.9521,  0.9856,  1.0203,  1.0568,  1.0946,
   1.1352,  1.1775,  1.2225,  1.2703,  1.3218,  1.3771,  1.4378,  1.5044,  1.5790,  1.6643,  1.7627,  1.8807,  2.0310,  2.2387,  2.5818,  3.3667, //15

  -3.6182, -2.9694, -2.6777, -2.5029, -2.3760, -2.2757, -2.1928, -2.1217, -2.0590, -2.0034, -1.9524, -1.9061, -1.8635, -1.8237, -1.7866, -1.7517, -1.7187,
  -1.6872, -1.6578, -1.6292, -1.6020, -1.4806, -1.3773, -1.2862, -1.2043,  0.4428,  0.5465,  0.6604,  0.7883,  0.9363,  0.9691,  1.0033,  1.0390,  1.0768,
   1.1162,  1.1577,  1.2017,  1.2485,  1.2988,  1.3529,  1.4115,  1.4756,  1.5478,  1.6301,  1.7244,  1.8376,  1.9796,  2.1742,  2.4956,  3.2078, //20

  -3.5489, -2.9338, -2.6541, -2.4850, -2.3620, -2.2652, -2.1845, -2.1150, -2.0538, -1.9992, -1.9496, -1.9043, -1.8626, -1.8236, -1.7870, -1.7523, -1.7200,
  -1.6892, -1.6598, -1.6319, -1.6052, -1.4850, -1.3824, -1.2917, -1.2102,  0.4344,  0.5376,  0.6513,  0.7788,  0.9264,  0.9589,  0.9929,  1.0284,  1.0654,
   1.1042,  1.1451,  1.1885,  1.2352,  1.2846,  1.3380,  1.3957,  1.4600,  1.5303,  1.6107,  1.7031,  1.8140,  1.9527,  2.1397,  2.4475,  3.1150, //25

  -3.5090, -2.9118, -2.6402, -2.4762, -2.3565, -2.2604, -2.1810, -2.1130, -2.0532, -1.9992, -1.9501, -1.9054, -1.8639, -1.8250, -1.7888, -1.7544, -1.7221,
  -1.6913, -1.6620, -1.6340, -1.6074, -1.4881, -1.3863, -1.2960, -1.2147,  0.4294,  0.5324,  0.6457,  0.7725,  0.9192,  0.9515,  0.9853,  1.0204,  1.0575,
   1.0962,  1.1370,  1.1801,  1.2262,  1.2756,  1.3288,  1.3863,  1.4495,  1.5186,  1.5973,  1.6889,  1.7965,  1.9325,  2.1164,  2.4142,  3.0575, //30

  -3.4538, -2.8850, -2.6232, -2.4629, -2.3456, -2.2530, -2.1756, -2.1085, -2.0488, -1.9957, -1.9476, -1.9034, -1.8624, -1.8244, -1.7889, -1.7552, -1.7233,
  -1.6932, -1.6645, -1.6368, -1.6106, -1.4925, -1.3910, -1.3014, -1.2202,  0.4217,  0.5246,  0.6375,  0.7640,  0.9097,  0.9422,  0.9756,  1.0109,  1.0475,
   1.0859,  1.1262,  1.1690,  1.2142,  1.2631,  1.3153,  1.3718,  1.4342,  1.5024,  1.5796,  1.6692,  1.7762,  1.9088,  2.0896,  2.3777,  2.9899, //40

  -3.4129, -2.8678, -2.6112, -2.4535, -2.3384, -2.2473, -2.1709, -2.1051, -2.0467, -1.9940, -1.9466, -1.9029, -1.8623, -1.8244, -1.7890, -1.7556, -1.7240,
  -1.6940, -1.6653, -1.6381, -1.6117, -1.4942, -1.3932, -1.3037, -1.2230,  0.4188,  0.5217,  0.6346,  0.7603,  0.9061,  0.9383,  0.9716,  1.0066,  1.0428,
   1.0813,  1.1214,  1.1638,  1.2092,  1.2573,  1.3092,  1.3655,  1.4270,  1.4955,  1.5726,  1.6606,  1.7655,  1.8970,  2.0726,  2.3564,  2.9504, //50

  -3.3676, -2.8431, -2.5968, -2.4432, -2.3307, -2.2404, -2.1650, -2.1000, -2.0425, -1.9913, -1.9446, -1.9014, -1.8608, -1.8238, -1.7890, -1.7557, -1.7244,
  -1.6945, -1.6664, -1.6393, -1.6132, -1.4961, -1.3956, -1.3068, -1.2260,  0.4143,  0.5171,  0.6296,  0.7555,  0.9004,  0.9323,  0.9656,  1.0001,  1.0361,
   1.0743,  1.1143,  1.1568,  1.2020,  1.2501,  1.3016,  1.3576,  1.4188,  1.4862,  1.5618,  1.6495,  1.7530,  1.8808,  2.0535,  2.3300,  2.9151, //75

  -3.3462, -2.8336, -2.5892, -2.4383, -2.3265, -2.2375, -2.1634, -2.0988, -2.0421, -1.9910, -1.9444, -1.9012, -1.8611, -1.8239, -1.7892, -1.7565, -1.7254,
  -1.6959, -1.6676, -1.6405, -1.6143, -1.4978, -1.3975, -1.3088, -1.2284,  0.4114,  0.5144,  0.6268,  0.7524,  0.8970,  0.9289,  0.9621,  0.9966,  1.0326,
   1.0705,  1.1106,  1.1528,  1.1974,  1.2452,  1.2962,  1.3517,  1.4121,  1.4793,  1.5554,  1.6418,  1.7454,  1.8726,  2.0433,  2.3199,  2.8880, //100

  -3.3259, -2.8200, -2.5810, -2.4314, -2.3210, -2.2333, -2.1593, -2.0960, -2.0396, -1.9888, -1.9425, -1.9001, -1.8606, -1.8239, -1.7890, -1.7565, -1.7256,
  -1.6962, -1.6682, -1.6414, -1.6153, -1.4990, -1.3992, -1.3104, -1.2304,  0.4089,  0.5113,  0.6237,  0.7492,  0.8933,  0.9255,  0.9586,  0.9930,  1.0291,
   1.0671,  1.1067,  1.1486,  1.1930,  1.2407,  1.2918,  1.3470,  1.4075,  1.4745,  1.5493,  1.6354,  1.7375,  1.8630,  2.0324,  2.3030,  2.8613, //150

  -3.3187, -2.8136, -2.5753, -2.4281, -2.3189, -2.2319, -2.1588, -2.0952, -2.0390, -1.9883, -1.9423, -1.9000, -1.8605, -1.8233, -1.7890, -1.7563, -1.7256,
  -1.6962, -1.6684, -1.6416, -1.6158, -1.4998, -1.3999, -1.3113, -1.2313,  0.4073,  0.5095,  0.6218,  0.7477,  0.8919,  0.9236,  0.9569,  0.9914,  1.0273,
   1.0649,  1.1047,  1.1467,  1.1910,  1.2384,  1.2892,  1.3441,  1.4042,  1.4704,  1.5455,  1.6312,  1.7329,  1.8586,  2.0276,  2.2977,  2.8536,

  -3.3051, -2.8137, -2.5759, -2.4289, -2.3199, -2.2322, -2.1583, -2.0949, -2.0382, -1.9875, -1.9417, -1.8995, -1.8607, -1.8242, -1.7897, -1.7573, -1.7265,
  -1.6970, -1.6690, -1.6422, -1.6165, -1.5005, -1.4009, -1.3125, -1.2326,  0.4065,  0.5091,  0.6214,  0.7467,  0.8903,  0.9224,  0.9552,  0.9895,  1.0253,
   1.0627,  1.1025,  1.1443,  1.1887,  1.2363,  1.2870,  1.3418,  1.4019,  1.4681,  1.5423,  1.6276,  1.7291,  1.8553,  2.0237,  2.2924,  2.8395,

  -3.2997, -2.8061, -2.5701, -2.4239, -2.3157, -2.2291, -2.1562, -2.0930, -2.0371, -1.9870, -1.9411, -1.8991, -1.8600, -1.8233, -1.7890, -1.7565, -1.7258,
  -1.6966, -1.6687, -1.6416, -1.6160, -1.5005, -1.4008, -1.3127, -1.2328,  0.4054,  0.5082,  0.6202,  0.7454,  0.8897,  0.9215,  0.9546,  0.9889,  1.0247,
   1.0624,  1.1019,  1.1435,  1.1878,  1.2347,  1.2855,  1.3398,  1.3999,  1.4655,  1.5402,  1.6257,  1.7273,  1.8518,  2.0192,  2.2854,  2.8322,

  -3.2900, -2.8030, -2.5694, -2.4235, -2.3156, -2.2298, -2.1571, -2.0938, -2.0381, -1.9872, -1.9416, -1.8994, -1.8603, -1.8238, -1.7894, -1.7568, -1.7262,
  -1.6970, -1.6691, -1.6425, -1.6167, -1.5015, -1.4018, -1.3138, -1.2337,  0.4040,  0.5069,  0.6193,  0.7441,  0.8883,  0.9202,  0.9531,  0.9875,  1.0232,
   1.0609,  1.1006,  1.1426,  1.1868,  1.2346,  1.2856,  1.3405,  1.4000,  1.4665,  1.5405,  1.6259,  1.7268,  1.8514,  2.0181,  2.2840,  2.8320,

  -3.2864, -2.7999, -2.5662, -2.4229, -2.3155, -2.2289, -2.1562, -2.0936, -2.0376, -1.9873, -1.9416, -1.8994, -1.8604, -1.8241, -1.7898, -1.7575, -1.7268,
  -1.6977, -1.6697, -1.6430, -1.6175, -1.5022, -1.4025, -1.3144, -1.2343,  0.4048,  0.5074,  0.6198,  0.7448,  0.8889,  0.9208,  0.9535,  0.9878,  1.0239,
   1.0614,  1.1009,  1.1424,  1.1866,  1.2336,  1.2843,  1.3388,  1.3989,  1.4646,  1.5387,  1.6239,  1.7250,  1.8489,  2.0156,  2.2816,  2.8240};


const double DickeyFullerTest::CVTable_ARD[] = {
  -6.1728, -4.8563, -4.3427, -4.0541, -3.8538, -3.6994, -3.5751, -3.4691, -3.3782, -3.2976, -3.2260, -3.1609, -3.1022, -3.0484, -2.9978, -2.9510, -2.9071,
  -2.8657, -2.8268, -2.7900, -2.7552, -2.6001, -2.4713, -2.3608, -2.2633, -0.7119, -0.6179, -0.5134, -0.3960, -0.2595, -0.2292, -0.1974, -0.1644, -0.1299,
  -0.0930, -0.0546, -0.0139,  0.0303,  0.0768,  0.1279,  0.1836,  0.2433,  0.3110,  0.3881,  0.4772,  0.5851,  0.7205,  0.9090,  1.2323,  1.9914,  // 10

  -5.2336, -4.3409, -3.9642, -3.7444, -3.5882, -3.4673, -3.3673, -3.2817, -3.2073, -3.1422, -3.0831, -3.0287, -2.9793, -2.9341, -2.8917, -2.8519, -2.8148,
  -2.7791, -2.7453, -2.7132, -2.6830, -2.5474, -2.4334, -2.3345, -2.2459, -0.7656, -0.6733, -0.5711, -0.4560, -0.3226, -0.2929, -0.2618, -0.2295, -0.1954,
  -0.1593, -0.1218, -0.0820, -0.0395,  0.0057,  0.0547,  0.1079,  0.1666,  0.2323,  0.3054,  0.3906,  0.4916,  0.6203,  0.7972,  1.0849,  1.7257,  // 15

  -4.8765, -4.1350, -3.8106, -3.6175, -3.4789, -3.3702, -3.2803, -3.2034, -3.1356, -3.0750, -3.0215, -2.9724, -2.9270, -2.8848, -2.8452, -2.8083, -2.7736,
  -2.7405, -2.7093, -2.6793, -2.6507, -2.5243, -2.4173, -2.3231, -2.2388, -0.7911, -0.6996, -0.5989, -0.4851, -0.3515, -0.3221, -0.2915, -0.2594, -0.2259,
  -0.1907, -0.1531, -0.1139, -0.0723, -0.0276,  0.0211,  0.0730,  0.1306,  0.1943,  0.2669,  0.3512,  0.4514,  0.5753,  0.7454,  1.0240,  1.6287,  // 20

  -4.6989, -4.0165, -3.7208, -3.5416, -3.4119, -3.3118, -3.2278, -3.1547, -3.0912, -3.0349, -2.9846, -2.9376, -2.8946, -2.8546, -2.8171, -2.7819, -2.7482,
  -2.7167, -2.6865, -2.6577, -2.6307, -2.5087, -2.4056, -2.3144, -2.2327, -0.8056, -0.7154, -0.6153, -0.5026, -0.3705, -0.3411, -0.3107, -0.2787, -0.2454,
  -0.2101, -0.1731, -0.1338, -0.0922, -0.0477,  0.0003,  0.0523,  0.1094,  0.1730,  0.2453,  0.3274,  0.4254,  0.5484,  0.7128,  0.9842,  1.5729,  // 25

  -4.5788, -3.9518, -3.6686, -3.4995, -3.3770, -3.2792, -3.1977, -3.1285, -3.0677, -3.0129, -2.9634, -2.9182, -2.8765, -2.8379, -2.8015, -2.7675, -2.7354,
  -2.7048, -2.6756, -2.6477, -2.6209, -2.5020, -2.4004, -2.3110, -2.2310, -0.8151, -0.7256, -0.6264, -0.5132, -0.3814, -0.3519, -0.3219, -0.2901, -0.2571,
  -0.2222, -0.1855, -0.1469, -0.1056, -0.0613, -0.0135,  0.0381,  0.0949,  0.1579,  0.2288,  0.3112,  0.4100,  0.5302,  0.6966,  0.9632,  1.5302,  // 30

  -4.4417, -3.8676, -3.6060, -3.4472, -3.3307, -3.2382, -3.1603, -3.0944, -3.0365, -2.9845, -2.9376, -2.8941, -2.8544, -2.8163, -2.7812, -2.7484, -2.7171,
  -2.6878, -2.6594, -2.6327, -2.6069, -2.4910, -2.3929, -2.3065, -2.2278, -0.8288, -0.7399, -0.6411, -0.5286, -0.3981, -0.3693, -0.3391, -0.3077, -0.2745,
  -0.2395, -0.2032, -0.1646, -0.1232, -0.0789, -0.0317,  0.0195,  0.0757,  0.1386,  0.2094,  0.2900,  0.3865,  0.5064,  0.6694,  0.9330,  1.4861,  // 40

  -4.3679, -3.8235, -3.5693, -3.4173, -3.3048, -3.2155, -3.1405, -3.0758, -3.0194, -2.9683, -2.9219, -2.8790, -2.8400, -2.8036, -2.7694, -2.7369, -2.7064,
  -2.6776, -2.6503, -2.6243, -2.5990, -2.4857, -2.3887, -2.3031, -2.2260, -0.8364, -0.7476, -0.6488, -0.5372, -0.4063, -0.3773, -0.3473, -0.3156, -0.2825,
  -0.2477, -0.2109, -0.1722, -0.1307, -0.0865, -0.0390,  0.0119,  0.0685,  0.1309,  0.2019,  0.2832,  0.3788,  0.4972,  0.6595,  0.9215,  1.4687,  // 50

  -4.2748, -3.7604, -3.5208, -3.3747, -3.2665, -3.1812, -3.1099, -3.0479, -2.9931, -2.9443, -2.9002, -2.8592, -2.8215, -2.7860, -2.7528, -2.7217, -2.6925,
  -2.6644, -2.6377, -2.6120, -2.5875, -2.4774, -2.3830, -2.2992, -2.2236, -0.8450, -0.7573, -0.6588, -0.5477, -0.4181, -0.3893, -0.3593, -0.3283, -0.2951,
  -0.2602, -0.2240, -0.1854, -0.1451, -0.1011, -0.0544, -0.0032,  0.0532,  0.1146,  0.1844,  0.2649,  0.3600,  0.4805,  0.6393,  0.8977,  1.4358,  // 75

  -4.2214, -3.7298, -3.4980, -3.3565, -3.2505, -3.1667, -3.0962, -3.0355, -2.9816, -2.9334, -2.8897, -2.8497, -2.8123, -2.7776, -2.7450, -2.7141, -2.6846,
  -2.6571, -2.6305, -2.6054, -2.5811, -2.4721, -2.3782, -2.2956, -2.2205, -0.8493, -0.7611, -0.6633, -0.5517, -0.4222, -0.3936, -0.3634, -0.3319, -0.2997,
  -0.2652, -0.2287, -0.1900, -0.1487, -0.1047, -0.0575, -0.0063,  0.0498,  0.1111,  0.1809,  0.2608,  0.3551,  0.4735,  0.6323,  0.8853,  1.4156,  // 100

  -4.1862, -3.7032, -3.4774, -3.3387, -3.2348, -3.1526, -3.0834, -3.0241, -2.9716, -2.9242, -2.8812, -2.8416, -2.8049, -2.7705, -2.7382, -2.7077, -2.6791,
  -2.6519, -2.6256, -2.6008, -2.5769, -2.4694, -2.3770, -2.2952, -2.2209, -0.8553, -0.7679, -0.6700, -0.5592, -0.4306, -0.4022, -0.3720, -0.3407, -0.3078,
  -0.2734, -0.2372, -0.1986, -0.1576, -0.1141, -0.0671, -0.0164,  0.0392,  0.1012,  0.1707,  0.2508,  0.3451,  0.4634,  0.6217,  0.8745,  1.4054,  // 150

  -4.1587, -3.6859, -3.4616, -3.3261, -3.2248, -3.1442, -3.0760, -3.0177, -2.9651, -2.9186, -2.8760, -2.8371, -2.8006, -2.7668, -2.7351, -2.7047, -2.6761,
  -2.6489, -2.6229, -2.5982, -2.5744, -2.4671, -2.3750, -2.2935, -2.2194, -0.8562, -0.7686, -0.6711, -0.5604, -0.4311, -0.4026, -0.3727, -0.3417, -0.3091,
  -0.2748, -0.2384, -0.2003, -0.1599, -0.1163, -0.0698, -0.0191,  0.0361,  0.0972,  0.1669,  0.2475,  0.3430,  0.4603,  0.6176,  0.8710,  1.3933,  // 200

  -4.1353, -3.6751, -3.4561, -3.3198, -3.2183, -3.1371, -3.0700, -3.0116, -2.9595, -2.9135, -2.8708, -2.8322, -2.7961, -2.7625, -2.7308, -2.7010, -2.6729,
  -2.6461, -2.6207, -2.5959, -2.5722, -2.4658, -2.3743, -2.2929, -2.2193, -0.8597, -0.7726, -0.6748, -0.5639, -0.4350, -0.4061, -0.3758, -0.3446, -0.3122,
  -0.2776, -0.2413, -0.2030, -0.1622, -0.1188, -0.0721, -0.0210,  0.0351,  0.0967,  0.1661,  0.2454,  0.3405,  0.4578,  0.6164,  0.8683,  1.3894,  // 300

  -4.1168, -3.6611, -3.4423, -3.3084, -3.2097, -3.1314, -3.0651, -3.0073, -2.9562, -2.9099, -2.8680, -2.8293, -2.7931, -2.7599, -2.7286, -2.6988, -2.6708,
  -2.6440, -2.6183, -2.5937, -2.5702, -2.4647, -2.3736, -2.2928, -2.2195, -0.8623, -0.7747, -0.6769, -0.5663, -0.4371, -0.4080, -0.3779, -0.3468, -0.3141,
  -0.2797, -0.2434, -0.2050, -0.1647, -0.1213, -0.0744, -0.0236,  0.0319,  0.0930,  0.1624,  0.2417,  0.3363,  0.4534,  0.6111,  0.8627,  1.3812,  // 500

  -4.0982, -3.6539, -3.4378, -3.3035, -3.2051, -3.1266, -3.0606, -3.0035, -2.9529, -2.9066, -2.8648, -2.8258, -2.7902, -2.7569, -2.7259, -2.6963, -2.6686,
  -2.6420, -2.6166, -2.5923, -2.5686, -2.4625, -2.3718, -2.2912, -2.2182, -0.8630, -0.7754, -0.6780, -0.5669, -0.4380, -0.4090, -0.3793, -0.3482, -0.3158,
  -0.2813, -0.2452, -0.2068, -0.1662, -0.1226, -0.0757, -0.0250,  0.0302,  0.0906,  0.1606,  0.2396,  0.3339,  0.4504,  0.6087,  0.8628,  1.3790,  // 1000

  -4.0932, -3.6446, -3.4307, -3.2987, -3.2005, -3.1217, -3.0563, -2.9990, -2.9483, -2.9028, -2.8610, -2.8227, -2.7875, -2.7544, -2.7231, -2.6940, -2.6663,
  -2.6395, -2.6142, -2.5899, -2.5666, -2.4615, -2.3713, -2.2904, -2.2174, -0.8630, -0.7756, -0.6782, -0.5669, -0.4378, -0.4094, -0.3796, -0.3481, -0.3153,
  -0.2814, -0.2455, -0.2073, -0.1672, -0.1238, -0.0778, -0.0275,  0.0276,  0.0898,  0.1593,  0.2395,  0.3342,  0.4515,  0.6093,  0.8588,  1.3761   // 10000
};


//#ifdef DEBUG_MyDickeyFullerTest
#ifdef TEST
    //unit test class . . . not shared, just declare/implement here
    class UTDickeyFullerTest : public UnitTest {
    public:
        UTDickeyFullerTest() : UnitTest(UnitTest::AREA_MATH, "DickeyFullerTest") {}

        bool run() {
        	double eps_mult = 1e+8;
        	*m_status << "\n\nTesting DickeyFullerTest...\n" << endl;
        	*m_status << "\n\nFor equality tests, we're using MACHEPS: " <<
        	eps_mult * std::numeric_limits<double>::epsilon()
        	<< endl;

        	size_t num_tests=0, failed_tests=0, passed_tests=0, fail_flag=0, main_tests=0, passed_main_tests=0;
        	size_t my_type_errors=0, their_type_errors=0;

        	int CORRECT_RESULT=-1;
        	int nrows=120, ncols=1;  // number of rows/cols in the matlab test data files.
        	int nres = 7; // There are 7 results that we will test. They are:
        	const string RESULT_NAME[] = { "H", "p-value", "t-stat", "c-value", "coeff", "se", "MSE" };
        	// h, pValue, stat, cValue, reg.coeff(1), reg.se(1), reg.MSE
        	// They are stored in the files DickeyFullerTestResults00?.txt

        	string matlab_test_datafile, matlab_results_file;

        	double *y=new double[nrows];  // y stores the data vector that we will test for unit root.
        	double *their_result = new double[nres];
			double *my_result = new double[nres+2]; // added 2 to store CVTable sample size and sig level values

			DickeyFullerTest *dftest;

        	for (int i=0; i<20; i++){
            	if (i<10){
            		CORRECT_RESULT=0;
            		matlab_test_datafile = DATA_DIR + DIR_SEPARATOR + "testing"
            				 + DIR_SEPARATOR + "DickeyFullerTestData00" + itoa(i,10) + ".txt";
            		matlab_results_file = DATA_DIR + DIR_SEPARATOR + "testing"
            				 + DIR_SEPARATOR + "DickeyFullerTestResults00" + itoa(i,10) + ".txt";
            	} else {
            		CORRECT_RESULT=1;
            		matlab_test_datafile = DATA_DIR + DIR_SEPARATOR + "testing"
            				 + DIR_SEPARATOR + "DickeyFullerTestData0" + itoa(i,10) + ".txt";
            		matlab_results_file = DATA_DIR + DIR_SEPARATOR + "testing"
            				 + DIR_SEPARATOR + "DickeyFullerTestResults0" + itoa(i,10) + ".txt";
            	}
            	*m_status << "\n>>> Test Set " << i << ": " << matlab_test_datafile << endl;

            	int cnt = read_array(nrows, ncols, y, matlab_test_datafile);
        		int rcnt = read_array(1, nres, their_result, matlab_results_file);

    			dftest = new DickeyFullerTest(nrows,y);
    			//dftest.display_attributes(m_status);

    			dftest->run_test();

    			// h, pValue, stat, cValue, reg.coeff(1), reg.se(1), reg.MSE
    			my_result[0] = (double) dftest->get_test_H();
    			my_result[1] = dftest->get_p_value();
    			my_result[2] = dftest->get_test_statistic();
    			my_result[3] = dftest->get_c_value();
    			my_result[4] = dftest->get_reg_coef((size_t)0);
    			my_result[5] = dftest->get_reg_se(0);
    			my_result[6] = dftest->get_reg_mse();
    			my_result[7] = (double) dftest->get_CV_sample_size();
    			my_result[8] = dftest->get_sig_level();

    			*m_status << "\t CV table indices: (sample size, sig level) = (" << my_result[7] << ", " << my_result[8] << ")\n" << endl;

    			for (int j=0; j<nres; j++){
            		num_tests++;
					if (j==0) {
						main_tests++;
						*m_status << "MAIN TEST:" << endl;
					}
    				if (equal(my_result[j],their_result[j],eps_mult)){
    					// test passed
    					if (j==0) {
    						*m_status << "\t\t ...PASSED: "  << RESULT_NAME[j] << "\t my_result[" << j << "] = " << my_result[j];
    						*m_status << " == " << their_result[j] << " = their_result[" << j << "]" << endl;
    						passed_main_tests++;
    					}
    					passed_tests++;
    				} else {
    					*m_status << "\t\t ...FAILED: " << RESULT_NAME[j] << "\t my_result[" << j << "] = " << my_result[j];
    					*m_status << " != \t" << their_result[j] << " = their_result[" << j << "]" << endl;
    					failed_tests++;
    				}
					if (j==0) {
						if (my_result[j]!=CORRECT_RESULT) {
							*m_status << "  WARNING: we've made a TYPE I or II error..." << endl;
							my_type_errors++;
						}
						if (their_result[j]!=CORRECT_RESULT) {
							*m_status << "\t\t\t\t\t\t ...matlab made a TYPE I or II error.\n" << endl;
							their_type_errors++;
						}
						*m_status << "OTHERS (failed):" << endl;
					}
    			}
    			delete dftest;
        	}
        	*m_status << "\nThe sample size assumed when estimating critical values: " << my_result[7] << endl;
        	*m_status << "\nDickeyFullerTest TEST RESULTS SUMMARY: "<< endl;
        	*m_status << "   PASSED " << passed_tests << "/" << num_tests << " TESTS" << endl;
        	*m_status << "   PASSED " << passed_main_tests << "/" << main_tests << " CRITICALLY IMPORTANT TESTS" << endl;
        	*m_status << "   TOTAL TYPE I or II ERRORS:  US " << my_type_errors << "/" << main_tests <<
						                         "     THEM " << their_type_errors << "/" << main_tests << endl;

    		delete[] y;	delete[] their_result; delete[] my_result;

    		if (passed_main_tests<main_tests){
        		return false;
            } else {
            	return true;
            }


        }
    };

    //to auto create/add to global test suite
    UTDickeyFullerTest dummy_UTDickeyFuller;
#endif

}
