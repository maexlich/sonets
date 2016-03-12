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
// Therefore, first calculate answer in that limit.
// Then use that estimate as an initial estimate for the 
// numerical routine to calculate for the general case.
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
  
  cout << "Before refine:\n";
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


  // status = calc_sqrtcov_given_rhos_refine
  //   (N, sigma, rho_recip, rho_conv, rho_div, rho_chain,
  //    sqrt_diag, sqrt_recip, sqrt_conv, sqrt_div, sqrt_chain);
  
  // cout << "After refine:\n";
  // cout << "sqrt_diag = ";
  // for(int i=0; i<2; i++)
  //   for(int j=0; j<2; j++)
  //     cout << sqrt_diag[i][j] << " ";
  // cout << "\n";
  // cout << "sqrt_recip = ";
  // for(int i=0; i<2; i++)
  //   for(int j=i; j<2; j++)
  //     cout << sqrt_recip[i][j] << " ";
  // cout << "\n";
  // cout << "sqrt_conv = ";
  // for(int i=0; i<2; i++)
  //   for(int j=0; j<2; j++)
  //     for(int k=j; k<2; k++)
  // 	cout << sqrt_conv[i][j][k] << " ";
  // cout << "\n";
  // cout << "sqrt_div = ";
  // for(int i=0; i<2; i++)
  //   for(int j=i; j<2; j++)
  //     for(int k=0; k<2; k++)
  // 	cout << sqrt_div[i][j][k] << " ";
  // cout << "\n";
  // cout << "sqrt_chain = ";
  // for(int i=0; i<2; i++)
  //   for(int j=0; j<2; j++)
  //     for(int k=0; k<2; k++)
  // 	cout << sqrt_chain[i][j][k] << " ";
  // cout << "\n";
  // cout.flush();

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
	  cout << "Found a negative eval(" << i <<")=" << the_eval 
	       << " for nrn_type=" << nrn_type << "\n";

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
      cout << "Can't calculate sqrt_diag[" << nrn_type << "][" << nrn_type 
	   << "]\n";
      gsl_eigen_symmv_free(work_eig);
      gsl_matrix_free(A);
      gsl_matrix_free(sqrtA);
      gsl_matrix_free(evecs);
      gsl_vector_free(evals);
      return -1;
    }
    
    sqrt_recip[nrn_type][nrn_type] = temp2/(2.0*sqrt_diag[nrn_type][nrn_type]);
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
      cout << "Couldn't calculate sqrt_diag[0][1] = " << sqrt_diag[0][1] <<"\n";
      return -1;
    }
    // find sign in front of temp3 to use for sqrt_diag[0][1]
    double temp_new = temp3*temp3/(temp12-sqrt_diag[0][1]*sqrt_diag[0][1]);
    if(temp_new < 0) {
      cout << "Couldn't calculate sqrt_diag[1][0] = " << sqrt_diag[1][0] <<"\n";
      return -1;
    }
    double thesign = GSL_SIGN((temp_new+temp12-temp21)
			      /(2*sqrt_diag[0][1]*sqrt(temp_new)));
    
    sqrt_diag[1][0] = thesign*sqrt(temp_new) - sqrt_diag[0][1];
    
    sqrt_recip[0][1] = temp3/(sqrt_diag[0][1]+sqrt_diag[1][0]);
  }


  return 0;

}


struct twopop_params
{
  double (*sigma)[2];
  double (*rho_recip)[2];
  double (*rho_conv)[2][2];
  double (*rho_div)[2][2];
  double (*rho_chain)[2][2];
  int *N;
};


int twopop_f(const gsl_vector * x, void *params,
	    gsl_vector * f) {
  
  double (*sigma)[2]= ((struct twopop_params *) params)->sigma;
  double (*rho_recip)[2]= ((struct twopop_params *) params)->rho_recip;
  double (*rho_conv)[2][2]= ((struct twopop_params *) params)->rho_conv;
  double (*rho_div)[2][2]= ((struct twopop_params *) params)->rho_div;
  double (*rho_chain)[2][2]= ((struct twopop_params *) params)->rho_chain;
  int *N= ((struct twopop_params *) params)->N;

  double a[2][2], b[2][2], c[2][2][2], d[2][2][2], e[2][2][2];
  
  int indcount=0;

  // four a's
  for(int i=0; i<2; i++) 
    for(int j=0; j<2; j++) {
      a[i][j] = gsl_vector_get(x,indcount);
      indcount++;
    }

  // three b's
  for(int i=0; i<2; i++) 
    for(int j=i; j<2; j++) {
      b[i][j] = gsl_vector_get(x,indcount);
      indcount++;
    }
  
  // six c's 
  for(int i=0; i<2; i++) 
    for(int j=0; j<2; j++)
      for(int k=j; k<2; k++) {
	c[i][j][k] = gsl_vector_get(x,indcount);
	indcount++;
      }

  // six d's 
  for(int i=0; i<2; i++) 
    for(int j=i; j<2; j++)
      for(int k=0; k<2; k++) {
	d[i][j][k] = gsl_vector_get(x,indcount);
	indcount++;
      }

  // eight e's
  for(int i=0; i<2; i++) 
    for(int j=0; j<2; j++)
      for(int k=0; k<2; k++) {
	e[i][j][k] = gsl_vector_get(x,indcount);
	indcount++;
      }
  
  double f_diag[2][2], f_recip[2][2], f_conv[2][2][2], f_div[2][2][2],
    f_chain[2][2][2];

  f_diag[0][0] =  gsl_pow_2(a[0][0]) + gsl_pow_2(b[0][0]) 
    + (N[0]-2)*(gsl_pow_2(c[0][0][0])+gsl_pow_2(d[0][0][0]) 
		+ 2*gsl_pow_2(e[0][0][0]))
    +N[1]*(gsl_pow_2(c[0][0][1])+gsl_pow_2(d[0][1][0]) 
	   + gsl_pow_2(e[0][0][1]) + gsl_pow_2(e[1][0][0]))
    - gsl_pow_2(sigma[0][0]);
  f_diag[0][1] = gsl_pow_2(a[0][1]) + gsl_pow_2(b[0][1]) 
    +(N[0]-1)*(gsl_pow_2(c[0][0][1])+gsl_pow_2(d[0][0][1]) 
	       + gsl_pow_2(e[0][1][0]) + gsl_pow_2(e[0][0][1]))
    +(N[1]-1)*(gsl_pow_2(c[0][1][1])+gsl_pow_2(d[0][1][1])
	       + gsl_pow_2(e[0][1][1]) + gsl_pow_2(e[1][0][1]))
    -gsl_pow_2(sigma[0][1]);
  f_diag[1][0] = gsl_pow_2(a[1][0]) + gsl_pow_2(b[0][1]) 
    +(N[0]-1)*(gsl_pow_2(c[1][0][0])+gsl_pow_2(d[0][1][0])
	       + gsl_pow_2(e[1][0][0]) + gsl_pow_2(e[0][1][0]))
    +(N[1]-1)*(gsl_pow_2(c[1][0][1])+gsl_pow_2(d[1][1][0])
	       + gsl_pow_2(e[1][0][1]) + gsl_pow_2(e[1][1][0]))
    - gsl_pow_2(sigma[1][0]);
  f_diag[1][1] =  gsl_pow_2(a[1][1]) + gsl_pow_2(b[1][1]) 
    +N[0]*(gsl_pow_2(c[1][0][1])+gsl_pow_2(d[0][1][1])
	   + gsl_pow_2(e[1][1][0]) + gsl_pow_2(e[0][1][1]))
    +(N[1]-2)*(gsl_pow_2(c[1][1][1])+gsl_pow_2(d[1][1][1]) 
	       + 2*gsl_pow_2(e[1][1][1]))
    - gsl_pow_2(sigma[1][1]);


  f_recip[0][0] =  2*a[0][0]*b[0][0]
    +(N[0]-2)*(2*c[0][0][0]*e[0][0][0] + 2*d[0][0][0]*e[0][0][0])
    +N[1]*(2*c[0][0][1]*e[0][0][1] + 2*d[0][1][0]*e[1][0][0] )
    -gsl_pow_2(sigma[0][0])*rho_recip[0][0];
  f_recip[0][1] =  (a[0][1]+a[1][0])*b[0][1]
    + (N[0]-1)*(c[0][0][1]*e[1][0][0] + c[1][0][0]*e[0][1][0]
		+ d[0][0][1]*e[0][1][0] +
		d[0][1][0]*e[0][0][1])
    + (N[1]-1)*(c[0][1][1]*e[1][0][1] + c[1][0][1]*e[0][1][1]
		+ d[0][1][1]*e[1][1][0] + d[1][1][0]*e[1][0][1])
    -sigma[0][1]*sigma[1][0]*rho_recip[0][1];
  f_recip[1][1] = 2*a[1][1]*b[1][1]
    +N[0]*(2*c[1][0][1]*e[1][1][0]  + 2*d[0][1][1]*e[0][1][1])
    +(N[1]-2)*(2*c[1][1][1]*e[1][1][1]  + 2*d[1][1][1]*e[1][1][1])
    -gsl_pow_2(sigma[1][1])*rho_recip[1][1];

  f_conv[0][0][0] =  2*a[0][0]*c[0][0][0]
    + 2*b[0][0]*e[0][0][0] + 2*d[0][0][0]*e[0][0][0]
    +(N[0]-3)*(gsl_pow_2(c[0][0][0]) + gsl_pow_2(e[0][0][0]))
    +N[1]*(gsl_pow_2(c[0][0][1]) + gsl_pow_2(e[1][0][0]))
    -gsl_pow_2(sigma[0][0])*rho_conv[0][0][0];
  f_conv[0][0][1] = (a[0][0]+a[0][1])*c[0][0][1]
    + b[0][0]*e[0][0][1]+b[0][1]*e[1][0][0] + d[0][1][0]*e[0][1][0] 
    + d[0][0][1]*e[0][0][1]
    +(N[0]-2)*(c[0][0][0]*c[0][0][1] + e[0][0][0]*e[0][0][1])
    +(N[1]-1)*(c[0][0][1]*c[0][1][1] + e[1][0][0]*e[1][0][1])
    -sigma[0][0]*sigma[0][1]*rho_conv[0][0][1];
  f_conv[0][1][1] = 2*a[0][1]*c[0][1][1]
    + 2*b[0][1]*e[1][0][1]+ 2*d[0][1][1]*e[0][1][1] 
    +(N[0]-1)*(gsl_pow_2(c[0][0][1]) + gsl_pow_2(e[0][0][1]))
    +(N[1]-2)*(gsl_pow_2(c[0][1][1]) + gsl_pow_2(e[1][0][1])) 
    -gsl_pow_2(sigma[0][1])*rho_conv[0][1][1];
  f_conv[1][0][0] = 2*a[1][0]*c[1][0][0]
    + 2*b[0][1]*e[0][1][0] +2*d[0][1][0]*e[1][0][0] 
    +(N[0]-2)*(gsl_pow_2(c[1][0][0]) + gsl_pow_2(e[0][1][0]))
    +(N[1]-1)*(gsl_pow_2(c[1][0][1]) + gsl_pow_2(e[1][1][0]))
    -gsl_pow_2(sigma[1][0])*rho_conv[1][0][0];
  f_conv[1][0][1] = (a[1][0]+a[1][1])*c[1][0][1]
    + b[0][1]*e[0][1][1]+b[1][1]*e[1][1][0] + d[1][1][0]*e[1][1][0]
    + d[0][1][1]*e[1][0][1]
    +(N[0]-1)*(c[1][0][0]*c[1][0][1] + e[0][1][0]*e[0][1][1])
    +(N[1]-2)*(c[1][0][1]*c[1][1][1] + e[1][1][0]*e[1][1][1])
    -sigma[1][0]*sigma[1][1]*rho_conv[1][0][1];
  f_conv[1][1][1] =2*a[1][1]*c[1][1][1]
    + 2*b[1][1]*e[1][1][1] + 2*d[1][1][1]*e[1][1][1]
    +N[0]*(gsl_pow_2(c[1][0][1]) + gsl_pow_2(e[0][1][1]))
    +(N[1]-3)*(gsl_pow_2(c[1][1][1]) + gsl_pow_2(e[1][1][1]))
    - gsl_pow_2(sigma[1][1])*rho_conv[1][1][1];

  f_div[0][0][0] = 2*a[0][0]*d[0][0][0]
    + 2*b[0][0]*e[0][0][0]  + 2*c[0][0][0]*e[0][0][0] 
    +(N[0]-3)*(gsl_pow_2(d[0][0][0]) + gsl_pow_2(e[0][0][0]))
    +N[1]*(gsl_pow_2(d[0][1][0]) + gsl_pow_2(e[0][0][1]))
    -gsl_pow_2(sigma[0][0])*rho_div[0][0][0];
  f_div[0][1][0] =(a[0][0]+a[1][0])*d[0][1][0]
    + b[0][0]*e[1][0][0] + b[0][1]*e[0][0][1] + c[0][0][1]*e[0][1][0]
    + c[1][0][0]*e[1][0][0]
    +(N[0]-2)*(d[0][0][0]*d[0][1][0] + e[0][0][0]*e[1][0][0])
    +(N[1]-1)*(d[0][1][0]*d[1][1][0] + e[0][0][1]*e[1][0][1]) 
    -sigma[0][0]*sigma[1][0]*rho_div[0][1][0];
  f_div[1][1][0] = 2*a[1][0]*d[1][1][0]
    +2*b[0][1]*e[1][0][1] + 2*c[1][0][1]*e[1][1][0] 
    +(N[0]-1)*(gsl_pow_2(d[0][1][0]) + gsl_pow_2(e[1][0][0]))
    +(N[1]-2)*(gsl_pow_2(d[1][1][0]) + gsl_pow_2(e[1][0][1]))
    -gsl_pow_2(sigma[1][0])*rho_div[1][1][0];
  f_div[0][0][1] = 2*a[0][1]*d[0][0][1]
    + 2*b[0][1]*e[0][1][0] + 2*c[0][0][1]*e[0][0][1] 
    +(N[0]-2)*(gsl_pow_2(d[0][0][1])+ gsl_pow_2(e[0][1][0]))
    +(N[1]-1)*(gsl_pow_2(d[0][1][1]) + gsl_pow_2(e[0][1][1])) 
    - gsl_pow_2(sigma[0][1])*rho_div[0][0][1];
  f_div[0][1][1] = (a[0][1]+a[1][1])*d[0][1][1]
    + b[0][1]*e[1][1][0] + b[1][1]*e[0][1][1] + c[0][1][1]*e[0][1][1] 
    + c[1][0][1]*e[1][0][1]
    +(N[0]-1)*(d[0][0][1]*d[0][1][1] + e[0][1][0]*e[1][1][0])
    +(N[1]-2)*(d[0][1][1]*d[1][1][1] + e[0][1][1]*e[1][1][1])
    -sigma[0][1]*sigma[1][1]*rho_div[0][1][1];
  f_div[1][1][1] = 2*a[1][1]*d[1][1][1]
    + 2*b[1][1]*e[1][1][1]+ 2*c[1][1][1]*e[1][1][1]
    +N[0]*(gsl_pow_2(d[0][1][1]) + gsl_pow_2(e[1][1][0]))
    +(N[1]-3)*(gsl_pow_2(d[1][1][1]) + gsl_pow_2(e[1][1][1]))
    - gsl_pow_2(sigma[1][1])*rho_div[1][1][1];

  f_chain[0][0][0] = 2*a[0][0]*e[0][0][0]
    + b[0][0]*(c[0][0][0] +d[0][0][0]) + c[0][0][0]*d[0][0][0] 
    + gsl_pow_2(e[0][0][0])
    +(N[0]-3)*(c[0][0][0]*e[0][0][0] + d[0][0][0]*e[0][0][0])
    +N[1]*(c[0][0][1]*e[0][0][1] + d[0][1][0]*e[1][0][0])
    -gsl_pow_2(sigma[0][0])*rho_chain[0][0][0];
  f_chain[0][0][1] = (a[0][0]+a[0][1])*e[0][0][1]
    + b[0][0]*c[0][0][1] + b[0][1]*d[0][1][0] + c[0][0][1]*d[0][0][1] 
    + e[1][0][0]*e[0][1][0]
    +(N[0]-2)*(c[0][0][1]*e[0][0][0] + d[0][0][0]*e[0][0][1])
    +(N[1]-1)*(c[0][1][1]*e[0][0][1] + d[0][1][0]*e[1][0][1])
    -sigma[0][0]*sigma[0][1]*rho_chain[0][0][1];
  f_chain[0][1][0] = (a[0][1]+a[1][0])*e[0][1][0]
    + b[0][1]*(c[1][0][0] + d[0][0][1]) + c[0][0][1]*d[0][1][0]
    + e[0][0][1]*e[1][0][0]
    +(N[0]-2)*(c[1][0][0]*e[0][1][0] + d[0][0][1]*e[0][1][0])
    +(N[1]-1)*(c[1][0][1]*e[0][1][1] + d[0][1][1]*e[1][1][0])
    -sigma[0][1]*sigma[1][0]*rho_chain[0][1][0];
  f_chain[0][1][1] = (a[0][1]+a[1][1])*e[0][1][1]
    + b[0][1]*c[1][0][1] + b[1][1]*d[0][1][1] + c[0][1][1]*d[0][1][1] 
    + e[1][0][1]*e[1][1][0]
    +(N[0]-1)*(c[1][0][1]*e[0][1][0] + d[0][0][1]*e[0][1][1])
    +(N[1]-2)*(c[1][1][1]*e[0][1][1] + d[0][1][1]*e[1][1][1])
    - sigma[0][1]*sigma[1][1]*rho_chain[0][1][1];
  f_chain[1][0][0] = (a[1][0]+a[0][0])*e[1][0][0]
    + b[0][1]*c[0][0][1] + b[0][0]*d[0][1][0] + c[1][0][0]*d[0][1][0] 
    + e[0][1][0]*e[0][0][1]
    +(N[0]-2)*(c[0][0][0]*e[1][0][0] + d[0][1][0]*e[0][0][0])
    +(N[1]-1)*(c[0][0][1]*e[1][0][1] + d[1][1][0]*e[1][0][0])
    - sigma[1][0]*sigma[0][0]*rho_chain[1][0][0];
  f_chain[1][0][1] = (a[1][0]+a[0][1])*e[1][0][1]
    + b[0][1]*(c[0][1][1] + d[1][1][0]) + c[1][0][1]*d[0][1][1] 
    + e[1][1][0]*e[0][1][1]
    +(N[0]-1)*(c[0][0][1]*e[1][0][0] + d[0][1][0]*e[0][0][1])
    +(N[1]-2)*(c[0][1][1]*e[1][0][1] + d[1][1][0]*e[1][0][1])
    - sigma[1][0]*sigma[0][1]*rho_chain[1][0][1];
  f_chain[1][1][0] =(a[1][1]+a[1][0])*e[1][1][0]
    + b[1][1]*c[1][0][1] + b[0][1]*d[0][1][1] + c[1][0][1]*d[1][1][0] 
    + e[0][1][1]*e[1][0][1]
    +(N[0]-1)*(c[1][0][0]*e[1][1][0] + d[0][1][1]*e[0][1][0])
    +(N[1]-2)*(c[1][0][1]*e[1][1][1] + d[1][1][1]*e[1][1][0]) 
    - sigma[1][1]*sigma[1][0]*rho_chain[1][1][0];
  f_chain[1][1][1] =  2*a[1][1]*e[1][1][1]
    + b[1][1]*(c[1][1][1] +d[1][1][1]) + c[1][1][1]*d[1][1][1] 
    + gsl_pow_2(e[1][1][1])
    +N[0]*(c[1][0][1]*e[1][1][0] + d[0][1][1]*e[0][1][1])
    +(N[1]-3)*(c[1][1][1]*e[1][1][1] + d[1][1][1]*e[1][1][1])
    - gsl_pow_2(sigma[1][1])*rho_chain[1][1][1];
  

  indcount=0;

  // four f_diag's
  for(int i=0; i<2; i++) 
    for(int j=0; j<2; j++) {
      gsl_vector_set(f,indcount, f_diag[i][j]);
      indcount++;
    }
  
  // three f_recip's
  for(int i=0; i<2; i++) 
    for(int j=i; j<2; j++) {
      gsl_vector_set(f,indcount, f_recip[i][j]);
      indcount++;
    }
  
  // six f_conv's
  for(int i=0; i<2; i++) 
    for(int j=0; j<2; j++)
      for(int k=j; k<2; k++) {
	gsl_vector_set(f,indcount, f_conv[i][j][k]);
	indcount++;
      }
  
  // six f_div's
  for(int i=0; i<2; i++) 
    for(int j=i; j<2; j++)
      for(int k=0; k<2; k++) {
	gsl_vector_set(f,indcount, f_div[i][j][k]);
	indcount++;
      }
  
  // eight f_chain's
  for(int i=0; i<2; i++) 
    for(int j=0; j<2; j++)
      for(int k=0; k<2; k++) {
	gsl_vector_set(f,indcount, f_chain[i][j][k]);
	indcount++;
      }

  return GSL_SUCCESS;
}


int calc_sqrtcov_given_rhos_refine
(int N[], double (*sigma)[2], double (*rho_recip)[2], double (*rho_conv)[2][2], 
 double (*rho_div)[2][2], double (*rho_chain)[2][2], 
 double (*sqrt_diag)[2], double (*sqrt_recip)[2], double (*sqrt_conv)[2][2], 
 double (*sqrt_div)[2][2], double (*sqrt_chain)[2][2]) {

  const gsl_multiroot_fsolver_type *T;
  gsl_multiroot_fsolver *s;
  
  int status;
  
  const size_t n = 27;
  
  struct twopop_params p = { sigma,rho_recip, 
  			    rho_conv, rho_div,
  			     rho_chain, N};
  

  gsl_multiroot_function f = {&twopop_f, n, &p};


  // set initial state based on the sqrt parameters passed into the function
  gsl_vector *x = gsl_vector_alloc (n);

  int indcount=0;
  // four sqrt_diag's
  for(int i=0; i<2; i++) 
    for(int j=0; j<2; j++) {
      gsl_vector_set(x,indcount, sqrt_diag[i][j]);
      indcount++;
    }

  // three sqrt_recip's
  for(int i=0; i<2; i++) 
    for(int j=i; j<2; j++) {
      gsl_vector_set(x,indcount, sqrt_recip[i][j]);
      indcount++;
    }
  
  // six sqrt_conv's
  for(int i=0; i<2; i++) 
    for(int j=0; j<2; j++)
      for(int k=j; k<2; k++) {
	gsl_vector_set(x,indcount, sqrt_conv[i][j][k]);
	indcount++;
      }

  // six sqrt_div's
  for(int i=0; i<2; i++) 
    for(int j=i; j<2; j++)
      for(int k=0; k<2; k++) {
	gsl_vector_set(x,indcount, sqrt_div[i][j][k]);
	indcount++;
      }

  // eight sqrt_chain's
  for(int i=0; i<2; i++) 
    for(int j=0; j<2; j++)
      for(int k=0; k<2; k++) {
	gsl_vector_set(x,indcount, sqrt_chain[i][j][k]);
	indcount++;
      }

  
  T = gsl_multiroot_fsolver_hybrids;
  s = gsl_multiroot_fsolver_alloc (T, n);
  gsl_multiroot_fsolver_set (s, &f, x);
  
  //print_state (iter, s, n);

  int iter=0;
  do
    {
      iter++;
      status = gsl_multiroot_fsolver_iterate (s);
      
      // print_state (iter, s, n);
      
      if (status)   /* check if solver is stuck */
	break;
      
      status =
	gsl_multiroot_test_residual (s->f, 1e-7);
    }
  while (status == GSL_CONTINUE && iter < 1000);
  
  // cout << "After " << iter << " iterations of solver, square root refine status = "
  //      << gsl_strerror(status) << "\n";

  // if had problems, don't modify values of sqrt parameter
  if(status) {
    gsl_multiroot_fsolver_free (s);
    gsl_vector_free (x);
    return(status);
  }  

  // else set sqrt parameters to result of function
  indcount=0;

  // four sqrt_diag's
  for(int i=0; i<2; i++) 
    for(int j=0; j<2; j++) {
      sqrt_diag[i][j] = gsl_vector_get(s->x,indcount);
      indcount++;
    }

  // three sqrt_recip's
  for(int i=0; i<2; i++) 
    for(int j=i; j<2; j++) {
      sqrt_recip[i][j] = gsl_vector_get(s->x,indcount);
      indcount++;
    }
  
  // six sqrt_conv's
  for(int i=0; i<2; i++) 
    for(int j=0; j<2; j++)
      for(int k=j; k<2; k++) {
	sqrt_conv[i][j][k] = gsl_vector_get(s->x,indcount);
	indcount++;
      }

  // six sqrt_div's
  for(int i=0; i<2; i++) 
    for(int j=i; j<2; j++)
      for(int k=0; k<2; k++) {
	sqrt_div[i][j][k] = gsl_vector_get(s->x,indcount);
	indcount++;
      }

  // eight sqrt_chain's
  for(int i=0; i<2; i++) 
    for(int j=0; j<2; j++)
      for(int k=0; k<2; k++) {
	sqrt_chain[i][j][k] = gsl_vector_get(s->x,indcount);
	indcount++;
      }


  gsl_multiroot_fsolver_free (s);
  gsl_vector_free (x);

  return 0;
}
