#include <stdlib.h>
#include <stdio.h>
#include <iostream>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_multiroots.h>
#include<gsl/gsl_eigen.h>
#include<gsl/gsl_blas.h>
#include<gsl/gsl_cdf.h>
#include "calc_sqrtcov_rec_2p.hpp"

using namespace std;

int calc_sqrtcov_given_rhos_large_N
(int N[], double (*sigma)[2], double (*rho_recip)[2], double (*rho_conv)[2][2], 
 double (*rho_div)[2][2], double (*rho_chain)[2][2], 
 double (*sqrt_diag)[2], double (*sqrt_recip)[2], double (*sqrt_conv)[2][2], 
 double (*sqrt_div)[2][2], double (*sqrt_chain)[2][2]);
int calc_sqrtcov_given_rhos_refine
(int N[], double (*sigma)[2], double (*rho_recip)[2], double (*rho_conv)[2][2], 
 double (*rho_div)[2][2], double (*rho_chain)[2][2], 
 double (*sqrt_diag)[2], double (*sqrt_recip)[2], double (*sqrt_conv)[2][2], 
 double (*sqrt_div)[2][2], double (*sqrt_chain)[2][2]);



//////////////////////////////////////////////////////////////////
// calc_sqrtcov_given_rhos
//
// calculate parameters determining the square root of 
// the gaussian covariance matrix
// given the paramters rho determining the covariance structure
//
// Since not all combinations of rho are valid
// it is important to be able to test for valid parameters.
// This is simple to do in the large N limit, as we can reduce
// most of the calculation to taking the square root 
// of small matrices. 
// Therefore, calculate answer in that limit.
//////////////////////////////////////////////////////////////////
int calc_sqrtcov_given_rhos
(int N[], double (*p)[2], double (*rho_recip)[2], double (*rho_conv)[2][2], 
 double (*rho_div)[2][2], double (*rho_chain)[2][2], 
 double (*sqrt_diag)[2], double (*sqrt_recip)[2], double (*sqrt_conv)[2][2], 
 double (*sqrt_div)[2][2], double (*sqrt_chain)[2][2]) {

  double sigma[2][2];
  for(int i=0; i<2; i++)
    for(int j=0; j<2; j++) 
       sigma[i][j] = 1.0/gsl_cdf_ugaussian_Qinv(p[i][j]);

  int status = calc_sqrtcov_given_rhos_large_N
    (N, sigma, rho_recip, rho_conv, rho_div, rho_chain,
     sqrt_diag, sqrt_recip, sqrt_conv, sqrt_div, sqrt_chain);
  
  if(status)
    return status;
  
  cout << "sqrt_diag = ";
  for(int i=0; i<2; i++)
    for(int j=0; j<2; j++)
      cout << sqrt_diag[i][j] << " ";
  cout << "\n";
  cout << "sqrt_recip = ";
  for(int i=0; i<2; i++)
    for(int j=i; j<2; j++)
      cout << sqrt_recip[i][j] << " ";
  cout << "\n";
  cout << "sqrt_conv = ";
  for(int i=0; i<2; i++)
    for(int j=0; j<2; j++)
      for(int k=j; k<2; k++)
  	cout << sqrt_conv[i][j][k] << " ";
  cout << "\n";
  cout << "sqrt_div = ";
  for(int i=0; i<2; i++)
    for(int j=i; j<2; j++)
      for(int k=0; k<2; k++)
  	cout << sqrt_div[i][j][k] << " ";
  cout << "\n";
  cout << "sqrt_chain = ";
  for(int i=0; i<2; i++)
    for(int j=0; j<2; j++)
      for(int k=0; k<2; k++)
  	cout << sqrt_chain[i][j][k] << " ";
  cout << "\n";
  cout.flush();

  // copy to make sqrt_recip symmetric
  for(int i=0; i<2; i++) 
    for(int j=i+1; j<2; j++) {
      sqrt_recip[j][i]=sqrt_recip[i][j];
    }
  
  // copy to make sqrt_conv symmetric in its last two components
  for(int i=0; i<2; i++) 
    for(int j=0; j<2; j++)
      for(int k=j+1; k<2; k++) {
	sqrt_conv[i][k][j]=sqrt_conv[i][j][k];
      }

  // copy to make sqrt_div symmetric in its first two components
  for(int i=0; i<2; i++) 
    for(int j=i+1; j<2; j++)
      for(int k=0; k<2; k++) {
	sqrt_div[j][i][k]=sqrt_div[i][j][k];
      }

  return status;

}


// calculate the components of the sqrt of the covariance matrix
// assuming that the number of neurons N[0] and N[1] are large
int calc_sqrtcov_given_rhos_large_N
(int N[], double (*sigma)[2], double (*rho_recip)[2], double (*rho_conv)[2][2], 
 double (*rho_div)[2][2], double (*rho_chain)[2][2], 
 double (*sqrt_diag)[2], double (*sqrt_recip)[2], double (*sqrt_conv)[2][2], 
 double (*sqrt_div)[2][2], double (*sqrt_chain)[2][2]) {

  int Ntot=N[0]+N[1];
  double eta[2] = {sqrt(N[0]/(double)Ntot), sqrt(N[1]/(double)Ntot)};

  // for large N, the equations for sqrt_conv, sqrt_div, and sqrt_chain
  // decouple from the rest.  More over, they decouple into two groups
  // of 10 equations based on the type of the central neuron in the motif
  
  // moreover, these 10 equations can be written as solving for the
  // square root of a 4x4 covariance matrix
  // Hence one can quickly determine if the equations have a real solution
  // and find that real solution
  
  const size_t nmat=4;
  gsl_eigen_symmv_workspace *work_eig=gsl_eigen_symmv_alloc(nmat);
  
  gsl_matrix *A = gsl_matrix_alloc(nmat,nmat);    // the 4x4 covariance matrix
  gsl_matrix *sqrtA =gsl_matrix_alloc(nmat,nmat); // its square root
  gsl_vector *evals=gsl_vector_alloc(nmat);       // evals of A
  gsl_matrix *evecs=gsl_matrix_alloc(nmat,nmat);  // evects of A

  for(int nrn_type=0; nrn_type <2; nrn_type++) {
    // Need to scale the rho's appropriately so the 10 equations
    // can be written as finding the square root of covariance matrix A
    gsl_matrix_set(A,0,0,sigma[nrn_type][0]*sigma[nrn_type][0]
		   *rho_conv[nrn_type][0][0]/(eta[1]*eta[1]));
    gsl_matrix_set(A,1,0,sigma[nrn_type][0]*sigma[nrn_type][1]
		   *rho_conv[nrn_type][0][1]/(eta[0]*eta[1]));
    gsl_matrix_set(A,2,0,sigma[0][nrn_type]*sigma[nrn_type][0]
		   *rho_chain[0][nrn_type][0]/(eta[1]*eta[1]));
    gsl_matrix_set(A,3,0,sigma[1][nrn_type]*sigma[nrn_type][0]
		   *rho_chain[1][nrn_type][0]/(eta[0]*eta[1]));
    gsl_matrix_set(A,1,1,sigma[nrn_type][1]*sigma[nrn_type][1]
		   *rho_conv[nrn_type][1][1]/(eta[0]*eta[0]));
    gsl_matrix_set(A,2,1,sigma[0][nrn_type]*sigma[nrn_type][1]
		   *rho_chain[0][nrn_type][1]/(eta[0]*eta[1]));
    gsl_matrix_set(A,3,1,sigma[1][nrn_type]*sigma[nrn_type][1]
		   *rho_chain[1][nrn_type][1]/(eta[0]*eta[0]));
    gsl_matrix_set(A,2,2,sigma[0][nrn_type]*sigma[0][nrn_type]
		   *rho_div[0][0][nrn_type]/(eta[1]*eta[1]));
    gsl_matrix_set(A,3,2,sigma[0][nrn_type]*sigma[1][nrn_type]
		   *rho_div[0][1][nrn_type]/(eta[0]*eta[1]));
    gsl_matrix_set(A,3,3,sigma[1][nrn_type]*sigma[1][nrn_type]
		   *rho_div[1][1][nrn_type]/(eta[0]*eta[0]));


    // cout << "A\n";
    // for(int i=0; i<nmat; i++) {
    //   for(int j=0; j<nmat; j++) {
    // 	cout << gsl_matrix_get(A,i,j) << " ";
    //   }
    //   cout << "\n";
    // }
    
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
	  cerr << "Found a negative eval(" << i <<")=" << the_eval 
	       << " for nrn_type=" << nrn_type << "\n";
	  cerr << "Cannot generate a network with combination of rho_conv, rho_div, and rho_chain centered around population " << nrn_type << ".\n";

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
		 
    // cout << "sqrt(A)\n";
    // for(int i=0; i<nmat; i++) {
    //   for(int j=0; j<nmat; j++) {
    // 	cout << gsl_matrix_get(sqrtA,i,j) << " ";
    //   }
    //   cout << "\n";
    // }
  
    // undo scaling to get elements of the square root of 
    // original covariance matrix
    sqrt_conv[nrn_type][0][0] = gsl_matrix_get(sqrtA,0,0)*eta[1]/eta[0]
      /sqrt(Ntot);
    sqrt_conv[nrn_type][0][1] = gsl_matrix_get(sqrtA,1,0)/sqrt(Ntot);
    sqrt_chain[0][nrn_type][0] = gsl_matrix_get(sqrtA,2,0)*eta[1]/eta[0]
      /sqrt(Ntot);
    sqrt_chain[1][nrn_type][0] = gsl_matrix_get(sqrtA,3,0)/sqrt(Ntot);
    sqrt_conv[nrn_type][1][1] = gsl_matrix_get(sqrtA,1,1)*eta[0]/eta[1]
      /sqrt(Ntot);
    sqrt_chain[0][nrn_type][1] = gsl_matrix_get(sqrtA,2,1)/sqrt(Ntot);
    sqrt_chain[1][nrn_type][1] = gsl_matrix_get(sqrtA,3,1)*eta[0]/eta[1]
      /sqrt(Ntot);
    sqrt_div[0][0][nrn_type] = gsl_matrix_get(sqrtA,2,2)*eta[1]/eta[0]
      /sqrt(Ntot);
    sqrt_div[0][1][nrn_type] = gsl_matrix_get(sqrtA,3,2)/sqrt(Ntot);
    sqrt_div[1][1][nrn_type] = gsl_matrix_get(sqrtA,3,3)*eta[0]/eta[1]
      /sqrt(Ntot);
    
    // can solve for sqrt_recip and sqrt_diag entries just involing nrn_type
    int onrn_type = 1-nrn_type;
    double temp1 = gsl_pow_2(sigma[nrn_type][nrn_type])
      - (N[nrn_type]-2)
      *(gsl_pow_2(sqrt_conv[nrn_type][nrn_type][nrn_type])
	+ gsl_pow_2(sqrt_div[nrn_type][nrn_type][nrn_type])
	+ 2*gsl_pow_2(sqrt_chain[nrn_type][nrn_type][nrn_type]))
      - N[onrn_type]*(gsl_pow_2(sqrt_conv[nrn_type][0][1])
		      + gsl_pow_2(sqrt_div[0][1][nrn_type])
		      + gsl_pow_2(sqrt_chain[0][nrn_type][1])
		      + gsl_pow_2(sqrt_chain[1][nrn_type][0]));

    // if temp1 is negative, will not be able to determine sqrt_diag
    // independent of value of rho_recip
    if(temp1 < 0) {
      cerr << "Can't calculate sqrt_diag[" << nrn_type << "][" << nrn_type 
	   << "]\n";
      cerr << "Cannot generate a network with combination of rho_conv, rho_div, and rho_chain  centered around population " << nrn_type << ".\n";
      
      gsl_eigen_symmv_free(work_eig);
      gsl_matrix_free(A);
      gsl_matrix_free(sqrtA);
      gsl_matrix_free(evecs);
      gsl_vector_free(evals);
      return -1;
    }


    double temp2 = gsl_pow_2(sigma[nrn_type][nrn_type])
      *rho_recip[nrn_type][nrn_type]
      - (N[nrn_type]-2)
      *(2*sqrt_conv[nrn_type][nrn_type][nrn_type]
	*sqrt_chain[nrn_type][nrn_type][nrn_type]
	+ 2*sqrt_div[nrn_type][nrn_type][nrn_type]
	*sqrt_chain[nrn_type][nrn_type][nrn_type])
      - N[onrn_type]*(2*sqrt_conv[nrn_type][0][1]
		      *sqrt_chain[nrn_type][nrn_type][onrn_type]
		      + 2*sqrt_div[0][1][nrn_type]
		      *sqrt_chain[onrn_type][nrn_type][nrn_type]);

    


    if(fabs(temp1) >= fabs(temp2)) {
      sqrt_diag[nrn_type][nrn_type] = sqrt((temp1+sqrt(temp1*temp1-temp2*temp2))
					   /2.0);
    }
    else {
      // if can't get real solution to sqrt_diag, original system did
      // not have a real solution (at least for large N)
      cerr << "Can't calculate sqrt_diag[" << nrn_type << "][" << nrn_type 
	   << "]\n";

      double sigma2 = gsl_pow_2(sigma[nrn_type][nrn_type]);
      double temp2a = rho_recip[nrn_type][nrn_type]*sigma2 - temp2;
      double rho_recip_max = GSL_MIN_DBL((fabs(temp1)+temp2a)/(sigma2),1);
      double rho_recip_min = GSL_MAX_DBL((-fabs(temp1)+temp2a)/(sigma2),-1);
    
      if(rho_recip_max > rho_recip_min) {
	cerr << "Cannot generate network when combine rho_recip with other parameters centered around population " << nrn_type << "\n";
	cerr << "Valid range of rho_recip given values of other parameters: ";
	cerr << "[" << rho_recip_min << ", " << rho_recip_max << "]\n";
	cerr << "(Specified value of rho_recip = " << rho_recip[nrn_type] << ")\n";
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
    
    if(sqrt_diag[nrn_type][nrn_type]) {
      sqrt_recip[nrn_type][nrn_type] = temp2/(2.0*sqrt_diag[nrn_type][nrn_type]);
    }
    else {
      sqrt_recip[nrn_type][nrn_type] = 0.0;
    }
  }

  gsl_eigen_symmv_free(work_eig);
  gsl_matrix_free(A);
  gsl_matrix_free(sqrtA);
  gsl_matrix_free(evecs);
  gsl_vector_free(evals);
  
//   cout << "c = ";
//   for(int i=0; i<2; i++)
//     for(int j=0; j<2; j++) 
//       for(int k=j; k<2; k++)
// 	cout << sqrt_conv[i][j][k] << ", ";
//   cout << "\n";
//   cout << "d = ";
//   for(int i=0; i<2; i++)
//     for(int j=0; j<2; j++) 
//       for(int k=j; k<2; k++)
// 	cout << sqrt_div[i][j][k] << ", ";
//   cout << "\n";
//   cout << "e = ";
//   for(int i=0; i<2; i++)
//     for(int j=0; j<2; j++) 
//       for(int k=0; k<2; k++)
// 	cout << sqrt_chain[i][j][k] << ", ";
//   cout << "\n";


  // lastly calculate sqrt_diag and sqrt_recip for terms involving
  // both neuron types
  double temp12 = gsl_pow_2(sigma[0][1])
    - (N[0]-1)*(gsl_pow_2(sqrt_conv[0][0][1])
		+ gsl_pow_2(sqrt_div[0][0][1])
		+ gsl_pow_2(sqrt_chain[0][1][0])
		+ gsl_pow_2(sqrt_chain[0][0][1]))
    - (N[1]-1)*(gsl_pow_2(sqrt_conv[0][1][1])
		+ gsl_pow_2(sqrt_div[0][1][1])
		+ gsl_pow_2(sqrt_chain[0][1][1])
		+ gsl_pow_2(sqrt_chain[1][0][1]));
  double temp21 = gsl_pow_2(sigma[1][0])
    - (N[0]-1)*(gsl_pow_2(sqrt_conv[1][0][0])
		+ gsl_pow_2(sqrt_div[0][1][0])
		+ gsl_pow_2(sqrt_chain[1][0][0])
		+ gsl_pow_2(sqrt_chain[0][1][0]))
    - (N[1]-1)*(gsl_pow_2(sqrt_conv[1][0][1])
		+ gsl_pow_2(sqrt_div[1][1][0])
		+ gsl_pow_2(sqrt_chain[1][0][1])
		+ gsl_pow_2(sqrt_chain[1][1][0]));
  double temp3 = sigma[0][1]*sigma[1][0]*rho_recip[0][1]
    - (N[0]-1)*(sqrt_conv[0][0][1]*sqrt_chain[1][0][0] 
		+ sqrt_conv[1][0][0]*sqrt_chain[0][1][0]
		+ sqrt_div[0][0][1]*sqrt_chain[0][1][0] 
		+ sqrt_div[0][1][0]*sqrt_chain[0][0][1])
    - (N[1]-1)*(sqrt_conv[0][1][1]*sqrt_chain[1][0][1] 
		+ sqrt_conv[1][0][1]*sqrt_chain[0][1][1]
		+ sqrt_div[0][1][1]*sqrt_chain[1][1][0] 
		+ sqrt_div[1][1][0]*sqrt_chain[1][0][1]);
  
  // cout << "temp12 = " << temp12 << "\n";
  // cout << "temp21 = " << temp21 << "\n";
  // cout << "temp3 = " << temp3 << "\n";

  if(fabs(temp3) < 1E-15) {
    sqrt_recip[0][1] = 0.0;
    sqrt_diag[0][1] = sqrt(temp12);
    sqrt_diag[1][0] = sqrt(temp21);
  }
  else {
    double temp4=temp12-temp21;
    double temp32=temp3*temp3;
    
    double tempa = temp4*temp4+4*temp32;
    double tempb = -4*temp12*temp32 - 2*temp32*temp4 - 2*temp12*temp4*temp4;
    double tempc = gsl_pow_2(temp32 + temp4*temp12);
    
    // cout << "tempa = " << tempa << ", tempb = " << tempb << ", tempc = " << tempc << endl;
    

    sqrt_diag[0][1] = sqrt((-tempb+sqrt(tempb*tempb-4*tempa*tempc))/(2*tempa));
    if(!(sqrt_diag[0][1] >= 0)) {
      cerr << "Couldn't calculate sqrt_diag[0][1] = " << sqrt_diag[0][1] <<"\n";
      return -1;
    }
    // find sign in front of temp3 to use for sqrt_diag[0][1]
    double temp_new = temp3*temp3/(temp12-sqrt_diag[0][1]*sqrt_diag[0][1]);
    if(temp_new < 0) {
      cerr << "Couldn't calculate sqrt_diag[1][0] = " << sqrt_diag[1][0] <<"\n";
      return -1;
    }
    double thesign = GSL_SIGN((temp_new+temp12-temp21)
			      /(2*sqrt_diag[0][1]*sqrt(temp_new)));
    
    sqrt_diag[1][0] = thesign*sqrt(temp_new) - sqrt_diag[0][1];

    if(sqrt_diag[0][1]+sqrt_diag[1][0]) {
      sqrt_recip[0][1] = temp3/(sqrt_diag[0][1]+sqrt_diag[1][0]);
    }
    else {
      sqrt_recip[0][1] = 0.0;
    }
  }


  return 0;

}
