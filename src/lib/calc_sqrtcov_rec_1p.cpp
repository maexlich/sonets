#include <stdlib.h>
#include <stdio.h>
#include <iostream>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_multiroots.h>
#include<gsl/gsl_eigen.h>
#include<gsl/gsl_blas.h>
#include<gsl/gsl_cdf.h>

#include "calc_sqrtcov_rec_1p.hpp"

using namespace std;

int calc_sqrtcov_given_rhos_large_N
(int N, double sigma, double rho_recip, double rho_conv, 
 double rho_div, double rho_chain,
 double &sqrt_diag, double &sqrt_recip, double &sqrt_conv,
 double &sqrt_div, double &sqrt_chain);


//////////////////////////////////////////////////////////////////
// calc_sqrtcov_given_rhos
//
// calculate parameters determining the square root of 
// the gaussian covariance matrix
// given the parameters rho determining the covariance structure
//
// Since not all combinations of rho are valid
// it is important to be able to test for valid parameters.
// This is simple to do in the large N limit, as we can reduce
// most of the calculation to taking the square root 
// of a small matrix. 
//////////////////////////////////////////////////////////////////
int calc_sqrtcov_given_rhos
(int N, double p, double rho_recip, double rho_conv, 
 double rho_div, double rho_chain, double rho_noshare,
 double &sqrt_diag, double &sqrt_recip, double &sqrt_conv,
 double &sqrt_div, double &sqrt_chain, double &sqrt_noshare) {

  // standard deviation of each component of the Gaussian
  // so that probability will be above zero will be p
  double sigma = 1.0/gsl_cdf_ugaussian_Qinv(p);

  // estimate standard deviations in the limit of a large network
  int status = calc_sqrtcov_given_rhos_large_N
    (N, sigma, rho_recip, rho_conv, rho_div, rho_chain,
     sqrt_diag, sqrt_recip, sqrt_conv, sqrt_div, sqrt_chain);
  sqrt_noshare=0.0;

  return status;

}


// calculate the components of the sqrt of the covariance matrix
// assuming that the number of neurons N is large
int calc_sqrtcov_given_rhos_large_N
(int N, double sigma, double rho_recip, double rho_conv, 
 double rho_div, double rho_chain,
 double &sqrt_diag, double &sqrt_recip, double &sqrt_conv,
 double &sqrt_div, double &sqrt_chain) {
  
  // for large N, the equations for sqrt_conv, sqrt_div, and sqrt_chain
  // decouple from the rest.
  
  // moreover, these 3 equations can be written as solving for the
  // square root of a 2x2 covariance matrix
  // Hence one can quickly determine if the equations have a real solution
  // and find that real solution

  const size_t nmat=2;
  gsl_eigen_symmv_workspace *work_eig=gsl_eigen_symmv_alloc(nmat);
  
  gsl_matrix *A = gsl_matrix_alloc(nmat,nmat);    // the 2x2 covariance matrix
  gsl_matrix *sqrtA =gsl_matrix_alloc(nmat,nmat); // its square root
  gsl_vector *evals=gsl_vector_alloc(nmat);       // evals of A
  gsl_matrix *evecs=gsl_matrix_alloc(nmat,nmat);  // evects of A

  gsl_matrix_set(A,0,0, rho_conv);
  gsl_matrix_set(A,1,1, rho_div);
  gsl_matrix_set(A,1,0, rho_chain);

  gsl_matrix_scale(A, sigma*sigma);

  // to calculate square root of A
  // 1. take it's eigen decomposition
  // 2. take the square root of its eigenvalues
  // 3. reconstruct with new eigenvalues to get a square root of A
  
  gsl_eigen_symmv(A, evals, evecs, work_eig);
  
  for(size_t i=0; i<nmat; i++) {
    double the_eval = gsl_vector_get(evals,i);
    if(the_eval <= 0) {
      if(the_eval > -1E-12) {
	// allow eigenvalues to be slightly negative due to
	// roundoff error
	the_eval=0;
      }
      else {
	// if have a negative eigenvalue, can't take square root
	// system of equations does not have a real solution
	// (at least in limit of large N)
	cerr << "Found a negative eval(" << i <<")=" << the_eval << "\n";
	cerr << "Cannot generate a network with combination of rho_conv, rho_div, and rho_chain.\n";
	gsl_eigen_symmv_free(work_eig);
	gsl_matrix_free(A);
	gsl_matrix_free(sqrtA);
	gsl_matrix_free(evecs);
	gsl_vector_free(evals);
	return -1;
      }
    }
    
    // scale eigenvector by fourth root of eval so 
    // reconstruction with be based on square root of eval
    gsl_vector_view col = gsl_matrix_column(evecs,i);
    gsl_vector_scale(&col.vector, sqrt(sqrt(the_eval)));
  }
  
  // reconstruct to get sqrt A
  gsl_blas_dgemm(CblasNoTrans,CblasTrans,1.0,evecs,evecs,0,sqrtA);
  
  // undo scaling to get elements of the square root of 
  // original covariance matrix
  // obtain approximate solutions for sqrt_conv, sqrt_div, and sqrt_chain
  gsl_matrix_scale(sqrtA,1.0/sqrt(N));
  sqrt_conv=gsl_matrix_get(sqrtA,0,0);
  sqrt_div=gsl_matrix_get(sqrtA,1,1);
  sqrt_chain=gsl_matrix_get(sqrtA,1,0);
  
  // can solve for sqrt_recip and sqrt_diag entries just involing nrn_type
  double temp1=sigma*sigma-(N-2.0)*(gsl_pow_2(sqrt_conv)
				    + gsl_pow_2(sqrt_div)
				    + 2.0*gsl_pow_2(sqrt_chain));

  // if temp1 is negative, will not be able to determine sqrt_diag
  // independent of value of rho_recip
  if(temp1 < 0) {
    cerr << "Can't calculate sqrt_diag\n";
    cerr << "Cannot generate a network with combination of rho_conv, rho_div, and rho_chain.\n";
    gsl_eigen_symmv_free(work_eig);
    gsl_matrix_free(A);
    gsl_matrix_free(sqrtA);
    gsl_matrix_free(evecs);
    gsl_vector_free(evals);
    return -1;
  }
   

  double temp2 = sigma*sigma*rho_recip
    -2*(N-2.0)*(sqrt_conv+sqrt_div)*sqrt_chain;
  
  // calculate sqrt_diag
  if(fabs(temp1) >= fabs(temp2)) {
    sqrt_diag = sqrt((temp1+sqrt(temp1*temp1-temp2*temp2))/2.0);
  }
  else {
    // if can't get real solution to sqrt_diag, original system did
    // not have a real solution (at least for large N)
    cerr << "Can't calculate sqrt_diag\n";

    double temp2a = 2*(N-2.0)*(sqrt_conv+sqrt_div)*sqrt_chain;
    double rho_recip_max = GSL_MIN_DBL((fabs(temp1)+temp2a)/(sigma*sigma),1);
    double rho_recip_min = GSL_MAX_DBL((-fabs(temp1)+temp2a)/(sigma*sigma),-1);
    
    if(rho_recip_max > rho_recip_min) {
      cerr << "Cannot generate network when combine rho_recip with other parameters\n";
      cerr << "Valid range of rho_recip given values of other parameters: ";
      cerr << "[" << rho_recip_min << ", " << rho_recip_max << "]\n";
      cerr << "(Specified value of rho_recip = " << rho_recip << ")\n";
    }
    else {
      cerr << "Cannot generate a network with combination of rho_conv, rho_div, and rho_chain.\n";
    }

    gsl_eigen_symmv_free(work_eig);
    gsl_matrix_free(A);
    gsl_matrix_free(sqrtA);
    gsl_matrix_free(evecs);
    gsl_vector_free(evals);
    return -1;
  }
  
  // calculate sqrt_recip
  sqrt_recip = temp2/(2.0*sqrt_diag);


  gsl_eigen_symmv_free(work_eig);
  gsl_matrix_free(A);
  gsl_matrix_free(sqrtA);
  gsl_matrix_free(evecs);
  gsl_vector_free(evals);

  return 0;

}
