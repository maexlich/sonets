#include <iostream>
#include <gsl/gsl_math.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_blas.h>

#include "calc_stats_2p.hpp"

using namespace std;


// calculate second order statistics of W containing two populations

int calc_N_sec(gsl_matrix *W, int N_nodes[2], 
	       double (*N_edge)[2],
	       double (*N_recip)[2], 
	       double (*N_conv)[2][2], double (*N_div)[2][2],
	       double (*N_chain)[2][2]) {

  

  // create the four submatrices
  gsl_matrix_view W00 = gsl_matrix_submatrix
    (W, 0, 0, N_nodes[0], N_nodes[0]);
  gsl_matrix_view W01 = gsl_matrix_submatrix
    (W, 0, N_nodes[0], N_nodes[0], N_nodes[1]);
  gsl_matrix_view W10 = gsl_matrix_submatrix
    (W, N_nodes[0], 0, N_nodes[1], N_nodes[0]);
  gsl_matrix_view W11 = gsl_matrix_submatrix
    (W, N_nodes[0], N_nodes[0], N_nodes[1], N_nodes[1]);

  // create a matrix of pointers to the submatrices
  gsl_matrix *Wsub[2][2];
  Wsub[0][0] = &W00.matrix;
  Wsub[0][1] = &W01.matrix;
  Wsub[1][0] = &W10.matrix;
  Wsub[1][1] = &W11.matrix;

  // calculate number of edge entries
  for(int ind1=0; ind1<2; ind1++)
    for(int ind2=0; ind2<2; ind2++) {
      N_edge[ind1][ind2] = 0.0;
      for(int k=0; k<N_nodes[ind1]; k++) {
	gsl_vector_view row = gsl_matrix_row(Wsub[ind1][ind2], k);
	N_edge[ind1][ind2] += gsl_blas_dasum(&row.vector);
      }
    }
  

  // create temporary matrices for to store multiplications of submatrices
  // need one of each size combination
  gsl_matrix *Wmult[2][2];
  for(int ind1=0; ind1<2; ind1++)
    for(int ind2=0; ind2<2; ind2++)
      Wmult[ind1][ind2] = gsl_matrix_alloc(N_nodes[ind1],N_nodes[ind2]);

  
  for(int ind1=0; ind1<2; ind1++)
    for(int ind2=0; ind2<2; ind2++) {

      // compute Wab*Wba, store in Wmult
      gsl_blas_dgemm (CblasNoTrans, CblasNoTrans, 1.0, 
		      Wsub[ind1][ind2], Wsub[ind2][ind1], 0.0, 
		      Wmult[ind1][ind1]);
      
      // first part of N_chain^aba is sum of Wab*Wba 
      N_chain[ind1][ind2][ind1] = 0.0;
      for(int i=0; i<N_nodes[ind1]; i++) {
	gsl_vector_view row = gsl_matrix_row(Wmult[ind1][ind1], i);
	N_chain[ind1][ind2][ind1] += gsl_blas_dasum(&row.vector);
      }

      if(ind2 >= ind1) {
	// N_recip is trace
	gsl_vector_view W2diag = gsl_matrix_diagonal(Wmult[ind1][ind1]);
	N_recip[ind1][ind2] = gsl_blas_dasum(&W2diag.vector);

	// need to subtract off N_chain
	N_chain[ind1][ind2][ind1] -= N_recip[ind1][ind2];
      }
      else {
	// this can only occur for ind1=2, ind2=1
	// so N_recip[1][2] has already been calculated
	// need to subtract it off N_chain
	N_chain[ind1][ind2][ind1] -= N_recip[ind2][ind1];
      }

      for(int ind3=0; ind3<2; ind3++) {
	
	// finish N_chain for ind3 != ind1
	if(ind3 != ind1) {
	  // compute Wab*Wbc, store in Wmult
	  gsl_blas_dgemm (CblasNoTrans, CblasNoTrans, 1.0, 
			  Wsub[ind1][ind2], Wsub[ind2][ind3], 0.0, 
			  Wmult[ind1][ind3]);
	  
	  // N_chain^abc is sum of Wab*Wbc
	  // don't have to subtract N_recip here
	  N_chain[ind1][ind2][ind3] = 0.0;
	  for(int i=0; i<N_nodes[ind1]; i++) {
	    gsl_vector_view row = gsl_matrix_row(Wmult[ind1][ind3], i);
	    N_chain[ind1][ind2][ind3] += gsl_blas_dasum(&row.vector);
	  }
	}

	if(ind3 >= ind2) {

	  // now use Wmult for Wab^T*Wac
	  gsl_blas_dgemm (CblasTrans, CblasNoTrans, 1.0, 
			  Wsub[ind1][ind2], Wsub[ind1][ind3], 0.0, 
			  Wmult[ind2][ind3]);
	  
	  // first part of N_conv^abc is sum of Wab^T*Wac 
	  N_conv[ind1][ind2][ind3]=0;
	  for(int i=0; i<N_nodes[ind2]; i++) {
	    gsl_vector_view row = gsl_matrix_row(Wmult[ind2][ind3], i);
	    N_conv[ind1][ind2][ind3] += gsl_blas_dasum(&row.vector);
	  }
	  
	  // have to subtract off N_edge if two converging pops are the same
	  if(ind2==ind3) 
	    N_conv[ind1][ind2][ind3] -= N_edge[ind1][ind2];
	  
	  // now use Wmult for Wba*Wca^T
	  gsl_blas_dgemm(CblasNoTrans, CblasTrans, 1.0, 
			 Wsub[ind2][ind1], Wsub[ind3][ind1], 0.0, 
			 Wmult[ind2][ind3]);
    
	  //first part of N_div^abc is sum of Wba*Wca^T
	  N_div[ind2][ind3][ind1] = 0.0;
	  for(int i=0; i<N_nodes[ind2]; i++) {
	    gsl_vector_view row = gsl_matrix_row(Wmult[ind2][ind3], i);
	    N_div[ind2][ind3][ind1] += gsl_blas_dasum(&row.vector);
	  }
	  
	  // have to subtract off N_edge if two diverging pops are the same
	  if(ind2==ind3) 
	    N_div[ind2][ind3][ind1] -= N_edge[ind2][ind1];
	}
      }
    }
  
    for(int ind1=0; ind1<2; ind1++)
      for(int ind2=0; ind2<2; ind2++)
	gsl_matrix_free(Wmult[ind1][ind2]);

  return 0;

}


int calc_phat_alphahat(gsl_matrix *W, int N_nodes[2], 
	       double (*phat)[2],
	       double (*alphahat_recip)[2], 
	       double (*alphahat_conv)[2][2], double (*alphahat_div)[2][2],
	       double (*alphahat_chain)[2][2]) {

  double N_edge[2][2];
  double N_recip[2][2];
  double N_conv[2][2][2];
  double N_div[2][2][2];
  double N_chain[2][2][2];
    
  calc_N_sec(W, N_nodes, N_edge, N_recip, N_conv, N_div, N_chain);

  
  for(int ind1=0; ind1<2; ind1++)
    for(int ind2=0; ind2<2; ind2++) {
      double denom = N_nodes[ind1]*N_nodes[ind2];
      if(ind1==ind2)
	denom -= N_nodes[ind1];
      phat[ind1][ind2] = N_edge[ind1][ind2]/denom;
    }

  for(int ind1=0; ind1<2; ind1++)
    for(int ind2=0; ind2<2; ind2++) {
      double denom = N_nodes[ind1]*N_nodes[ind2];
      if(ind1==ind2)
	denom -= N_nodes[ind1];
      if(ind2 >= ind1) 
	alphahat_recip[ind1][ind2] = N_recip[ind1][ind2]
	  / (denom*phat[ind1][ind2]*phat[ind2][ind1]) - 1.0;

      for(int ind3=0; ind3<2; ind3++) {
	denom = N_nodes[ind1]*(double) N_nodes[ind2]*(double)N_nodes[ind3];
	if(ind1 == ind2) {
	  if(ind1 == ind3) 
	    denom += -3*gsl_pow_2(N_nodes[ind1]) + 2*N_nodes[ind1];
	  else
	    denom -= N_nodes[ind1]*N_nodes[ind3];
	}
	else if(ind1 == ind3)
	  denom -= N_nodes[ind1]*N_nodes[ind2];
	else if(ind2 == ind3)
	  denom -= N_nodes[ind1]*N_nodes[ind2];

	if(ind3 >= ind2) {
	  alphahat_conv[ind1][ind2][ind3] = N_conv[ind1][ind2][ind3]
	    / (denom*phat[ind1][ind2]*phat[ind1][ind3]) - 1.0;
	  alphahat_div[ind2][ind3][ind1] = N_div[ind2][ind3][ind1]
	    / (denom*phat[ind2][ind1]*phat[ind3][ind1]) - 1.0;
	}

	alphahat_chain[ind1][ind2][ind3] = N_chain[ind1][ind2][ind3]
	  / (denom*phat[ind1][ind2]*phat[ind2][ind3]) - 1.0;

      }
    }

  return 0;

}


int calc_gaus_covs(gsl_matrix *W, int N_nodes[2], 
		   double (*sigma)[2],
		   double (*cov_recip)[2], 
		   double (*cov_conv)[2][2], double (*cov_div)[2][2],
		   double (*cov_chain)[2][2],
		   double &cov_other) {
  
  int N_shift[2];
  N_shift[0]=0;
  N_shift[1]=N_nodes[0];

  int N_nodes_tot= N_nodes[0]+N_nodes[1];

  // calc covariances of W (assume everything mean zero)
  
  cov_other=0.0;
    
  for(int i=0; i<2; i++)
    for(int j=0; j<2; j++) {
      sigma[i][j]=0.0;
      cov_recip[i][j]=0.0;
      for(int k=0; k<2; k++) {
	cov_conv[i][j][k]=0.0;
	cov_div[i][j][k]=0.0;
	cov_chain[i][j][k]=0.0;
      }
    }
    
  for(int ntype1=0; ntype1<2; ntype1++) 
    for(int i=N_shift[ntype1]; i<N_nodes[ntype1]+N_shift[ntype1]; i++) {
      for(int ntype2=0; ntype2<2; ntype2++) 
	for(int j=N_shift[ntype2]; j<N_nodes[ntype2]+N_shift[ntype2]; j++) {
	    
	  if(i==j)
	    continue;
	    
	  double w_ij = gsl_matrix_get(W,i,j);
	  sigma[ntype1][ntype2] += w_ij*w_ij;
	    
	  if(ntype2 < ntype1)
	    cov_recip[ntype2][ntype1] += w_ij * gsl_matrix_get(W,j,i);
	  else
	    cov_recip[ntype1][ntype2] += w_ij * gsl_matrix_get(W,j,i);
	    
	  for(int ntype3=0; ntype3<2; ntype3++) 
	    for(int k=N_shift[ntype3]; k<N_nodes[ntype3]+N_shift[ntype3]; k++) {
	      if(k==i || k==j)
		continue;
	      
	      if(ntype3 < ntype2) {
		cov_conv[ntype1][ntype3][ntype2]
		  += w_ij * gsl_matrix_get(W, i,k);
	      }
	      else {
		cov_conv[ntype1][ntype2][ntype3]
		  += w_ij * gsl_matrix_get(W, i,k);
	      }
	      
	      if(ntype3 < ntype1)
		cov_div[ntype3][ntype1][ntype2]
		  += w_ij * gsl_matrix_get(W, k,j);
	      else
		cov_div[ntype1][ntype3][ntype2]
		  += w_ij * gsl_matrix_get(W, k,j);


	      cov_chain[ntype1][ntype2][ntype3] 
		+= w_ij * gsl_matrix_get(W, j,k);
	      cov_chain[ntype3][ntype1][ntype2]
		+= w_ij * gsl_matrix_get(W, k,i);
		
	      // subsample edges that don't share a node
	      if(k+1 <N_nodes_tot && k+1 != i && k+1 !=j) 
		cov_other += w_ij * gsl_matrix_get(W, k, k+1);
	    
	    }
	    
	}
    } 
    
  for(int i=0; i<2; i++)
    for(int j=0; j<2; j++) {
      if(i==j) {
	sigma[i][i] /= N_nodes[i]*(N_nodes[i]-1.0);
	sigma[i][i] = sqrt(sigma[i][i]);
      }
      else {
	sigma[i][j] /= N_nodes[i]*N_nodes[j];
	sigma[i][j] = sqrt(sigma[i][j]);
      }
    }

  for(int i=0; i<2; i++)
    for(int j=i; j<2; j++) {
      if(i==j) {
	cov_recip[i][i] /= N_nodes[i]*(N_nodes[i]-1.0)*sigma[i][i]*sigma[i][i];
      }
      else 
	cov_recip[i][j] /= 2*N_nodes[i]*N_nodes[j]*sigma[i][j]*sigma[j][i];
    }
  
  for(int i=0; i<2; i++)
    for(int j=0; j<2; j++) 
      for(int k=j; k<2; k++) {
	if(j==k) {
	  if(i==j) {
	    cov_conv[i][i][i] /= N_nodes[i]*(N_nodes[i]-1.0)*(N_nodes[i]-2.0)
	      *sigma[i][i]*sigma[i][i];
	  }
	  else {
	    cov_conv[i][j][j] /= N_nodes[j]*(N_nodes[j]-1.0)*N_nodes[i]
	      *sigma[i][j]*sigma[i][j];
	  }
	}
	else {
	  if(i==j) {
	    cov_conv[i][i][k] /= 2*N_nodes[i]*(N_nodes[i]-1.0)*N_nodes[k]
	      *sigma[i][i]*sigma[i][k];
	  }
	  else if(i==k) {
	    cov_conv[i][j][i] /= 2*N_nodes[i]*(N_nodes[i]-1.0)*N_nodes[j]
	      *sigma[i][j]*sigma[i][i];
	  }
	  else {
	      cov_conv[i][j][k] /= 2*N_nodes[i]*N_nodes[j]*N_nodes[k]
		*sigma[i][j]*sigma[i][k];
	  }
	}
      }
  
	
  for(int i=0; i<2; i++)
    for(int j=i; j<2; j++) 
      for(int k=0; k<2; k++) {
	if(j==k) {
	  if(i==j) {
	    cov_div[i][i][i] /= N_nodes[i]*(N_nodes[i]-1.0)*(N_nodes[i]-2.0)
	      *sigma[i][i]*sigma[i][i];
	  }
	  else {
	    cov_div[i][j][j] /= 2*N_nodes[j]*(N_nodes[j]-1.0)*N_nodes[i]
	      *sigma[i][j]*sigma[j][j];
	  }
	}
	else {
	  if(i==j) {
	    cov_div[i][i][k] /= N_nodes[i]*(N_nodes[i]-1.0)*N_nodes[k]
	      *sigma[i][k]*sigma[i][k];
	  }
	  else if(i==k) {
	      cov_div[i][j][i] /= 2*N_nodes[i]*(N_nodes[i]-1.0)*N_nodes[j]
		*sigma[j][i]*sigma[i][i];
	    }
	  else {
	    cov_div[i][j][k] /= 2*N_nodes[i]*N_nodes[j]*N_nodes[k]
	      *sigma[i][k]*sigma[j][k];
	  }
	}
	
      }
  
  for(int i=0; i<2; i++)
    for(int j=0; j<2; j++) 
      for(int k=0; k<2; k++) {
	if(j==k) {
	  if(i==j) {
	    cov_chain[i][i][i] /= 2*N_nodes[i]*(N_nodes[i]-1.0)
	      *(N_nodes[i]-2.0)*sigma[i][i]*sigma[i][i];
	  }
	  else {
	    cov_chain[i][j][j] /= 2*N_nodes[j]*(N_nodes[j]-1.0)*N_nodes[i]
	      *sigma[i][j]*sigma[j][j]; 
	  }
	}
	else {
	  if(i==j) 
	    cov_chain[i][i][k] /= 2*N_nodes[i]*(N_nodes[i]-1.0)*N_nodes[k]
	      *sigma[i][i]*sigma[i][k];
	  else if(i==k)
	    cov_chain[i][j][i] /= 2*N_nodes[i]*(N_nodes[i]-1.0)*N_nodes[j]
	      *sigma[i][j]*sigma[j][i];
	  else 
	    cov_chain[i][j][k] /= 2*N_nodes[i]*N_nodes[j]*N_nodes[k]
	      *sigma[i][j]*sigma[j][k];
	  
	}
      }

  cov_other /= (N_nodes_tot-1.0)*(N_nodes_tot-2.0)*(N_nodes_tot-3.0);

  return 0;
}
