#include <iostream>
#include <string.h>
#include <fcntl.h>
#include <sys/stat.h>
#include <gsl/gsl_rng.h>
#include "secorder_rec_2p.hpp"
#include "calc_stats_2p.hpp"


using namespace std;



int main(int argc, char *argv[]) {
  
  // optional flags
  int deterministic_seed = 0; // if nonzero, use for random num gen seed

  int N_nodes[2] = {4000,1000};
  double p[2][2] = {{1E-10, 0.1},{0.1,0.1}};
  double alpha_recip[2][2] = {{0.0, 0.0},{99,00}};


  // form of large N covariance matrix centered around population x
  // conv_x00 conv_x01 chain_0x0 chain_1x0
  //          conv_x11 chain_0x1 chain_1x1
  //                     div_00x   div_01x
  //                               div_11x


  ///////////////////////////////////////////////
  // correlations centered around population zero
  // diagonal entries
  double alpha_conv_hom_00 = 0;  // zeros onto zero (1,1)
  double alpha_conv_hom_01 = 0.3;  // ones onto zero  (2,2)
  double alpha_div_hom_00  = 0;  // zero onto zeros (3,3)
  double alpha_div_hom_10  = 0.3;  // zero onto ones  (4,4)

  // cross entries for divergence and convergence
  // constrained by the diagonal entries
  double cc_conv_0_01 = 0;  // mix onto zero  (2,1)
  double cc_div_01_0  = 0;  // zero onto mix  (4,3)

  // remaining four cross entries are chains through population 0
  // constrained by diagonal entries
  double cc_chain_000 = 0;  // (3,1)
  double cc_chain_100 = 0;  // (4,1)
  double cc_chain_001 = 0;  // (3,2)
  double cc_chain_101 = 0;  // (4,2)

  ///////////////////////////////////////////////
  // correlations centered around population one
  // diagonal entries
  double alpha_conv_hom_10 = 1.0;  // zeros onto one (1,1)
  double alpha_conv_hom_11 = 0.3;  // ones onto one  (2,2)
  double alpha_div_hom_01  = 1.0;  // one onto zeros (3,3)
  double alpha_div_hom_11  = 1.0;  // one onto ones  (4,4)

  // cross entries for divergence and convergence
  // constrained by the diagonal entries
  double cc_conv_1_01 = 0.0;  // mix onto one  (2,1)
  double cc_div_01_1  = -0.9;  // one onto mix  (4,3)

  // remaining four cross entries are chains through population 1
  // constrained by diagonal entries
  double cc_chain_010 = -0.9;  // (3,1)
  double cc_chain_111 = 0.0;  // (4,1)
  double cc_chain_110 = 0.9;  // (3,2)
  double cc_chain_011 = 0.0;  // (4,2)

  // order of the indices for the 3D arrays
  // {{{000,001},{010,011}},{{100,101},{110,111}}}
  // Put 99 in entries that are ignored
  double alpha_conv_hom[2][2] = {{alpha_conv_hom_00,alpha_conv_hom_01},{alpha_conv_hom_10,alpha_conv_hom_11}};
  double cc_conv_mixed[2][2][2] = {{{99,cc_conv_0_01},{99,99}},{{99,cc_conv_1_01},{99,99}}};
  double alpha_div_hom[2][2] = {{alpha_div_hom_00,alpha_div_hom_01},{alpha_div_hom_10,alpha_div_hom_11}};
  double cc_div_mixed[2][2][2] = {{{99,99},{cc_div_01_0,cc_div_01_1}},{{99,99},{99,99}}};
  double cc_chain[2][2][2] = {{{cc_chain_000,cc_chain_001},{cc_chain_010,cc_chain_011}},{{cc_chain_100,cc_chain_101},{cc_chain_110,cc_chain_111}}};




  gsl_rng *r = gsl_rng_alloc(gsl_rng_mt19937);

  // set the seed
  // if deterministic_seed, use that for seed
  if(deterministic_seed)
    gsl_rng_set(r,deterministic_seed); 
  // else use time in seconds for the seed
  else
    gsl_rng_set(r, time(NULL));

  
  // generate second order matrix W
  gsl_matrix *W;
  W = secorder_rec_2p(N_nodes,p, alpha_recip, alpha_conv_hom, cc_conv_mixed,
		      alpha_div_hom, cc_div_mixed, cc_chain,
		      r);
  
  // if failed to generate a matrix, write error and quit
  if(!W) {
    cerr << "Failed to generate the matrix\n";
    return -1;
  }
 

  // matrix files will be stored in data directory
  // with filename determined by the six input parameters
  mkdir("data",0755);  // make data directory, if it doesn't exist
  char FNbase[300];

  sprintf(FNbase,"_%i_%i_%1.3f_%1.3f_%1.3f_%1.3f_%1.2f_%1.2f_%1.2f_%1.2f_%1.2f_%1.2f_%1.2f_%1.2f_%1.2f_%1.2f_%1.2f_%1.2f_%1.2f_%1.2f_%1.2f_%1.2f_%1.2f_%1.2f_%1.2f_%1.2f_%1.2f_%1.2f_%1.2f",
	  N_nodes[0], N_nodes[1],
	  p[0][0],p[0][1],p[1][0],p[1][1],
	  alpha_recip[0][0], alpha_recip[0][1], alpha_recip[1][1],
	  alpha_conv_hom[0][0], alpha_conv_hom[0][1],
	  alpha_conv_hom[1][0], alpha_conv_hom[1][1],
	  cc_conv_mixed[0][0][1], cc_conv_mixed[1][0][1], 
	  alpha_div_hom[0][0], alpha_div_hom[0][1], 
	  alpha_div_hom[1][0], alpha_div_hom[1][1],
	  cc_div_mixed[0][1][0], cc_div_mixed[0][1][1], 
	  cc_chain[0][0][0], cc_chain[0][0][1], 
	  cc_chain[0][1][0], cc_chain[0][1][1],
	  cc_chain[1][0][0], cc_chain[1][0][1],
	  cc_chain[1][1][0], cc_chain[1][1][1]);


  // write matrix to a file
  char FN[200];
  FILE *fhnd;
  strcpy(FN, "data/w");
  strcat(FN, FNbase);
  strcat(FN, ".dat");
  fhnd = fopen(FN, "w");
  if(fhnd==NULL) {
    cerr << "Couldn't open outfile file " << FN << "\n";
    exit(-1);
  }


  for(int i=0; i<N_nodes[0]+N_nodes[1];i++) {
    for(int j=0; j<N_nodes[0]+N_nodes[1]; j++) {
      fprintf(fhnd, "%i ", (int) gsl_matrix_get(W,i,j));
    }
    fprintf(fhnd,"\n");
  }
  fclose(fhnd);



  ////////////////////////////////////////////////////////////
  // Calculate the covariance structure of the matrix 
  // This should approximately agree with the input alphas
  ////////////////////////////////////////////////////////////

  cout << "Testing statistics of matrix...";
  cout.flush();


  double phat[2][2];
  double alphahat_recip[2][2];
  double alphahat_conv[2][2][2];
  double alphahat_div[2][2][2];
  double alphahat_chain[2][2][2];

  calc_phat_alphahat(W, N_nodes, phat, alphahat_recip, 
		     alphahat_conv, alphahat_div, alphahat_chain);

  cout << "done\n";

  strcpy(FN, "data/stats");
  strcat(FN, FNbase);
  strcat(FN, ".dat");
  fhnd = fopen(FN, "w");
  if(fhnd==NULL) {
    cerr << "Couldn't open outfile file " << FN << "\n";
    exit(-1);
  }

  
  fprintf(fhnd, "%i %i ", N_nodes[0], N_nodes[1]);
  
  cout << "Actual statistics of matrix:\n";
  cout << "phat = ";
  for(int i=0; i<2; i++)
    for(int j=0; j<2; j++) {
      cout << phat[i][j] << " ";
      fprintf(fhnd, "%e ", phat[i][j]);
    }
  cout << "\n";
  cout << "alphahat_recip = ";
  for(int i=0; i<2; i++)
    for(int j=i; j<2; j++) {
      cout << alphahat_recip[i][j] << " ";
      fprintf(fhnd, "%e ", alphahat_recip[i][j]);
    }
  cout << "\n";
  cout << "alphahat_conv = "; 
  for(int i=0; i<2; i++)
    for(int j=0; j<2; j++)
      for(int k=j; k<2; k++) {
	cout << alphahat_conv[i][j][k] << " ";
	fprintf(fhnd, "%e ", alphahat_conv[i][j][k]);
      }
  cout << "\n";
  cout << "alphahat_div = ";
  for(int i=0; i<2; i++)
    for(int j=i; j<2; j++)
      for(int k=0; k<2; k++) {
	cout << alphahat_div[i][j][k] << " ";
	fprintf(fhnd, "%e ", alphahat_div[i][j][k]);
      }
  cout << "\n";
  cout << "alphahat_chain = ";
  for(int i=0; i<2; i++)
    for(int j=0; j<2; j++)
      for(int k=0; k<2; k++) {
	cout << alphahat_chain[i][j][k] << " ";
	fprintf(fhnd, "%e ", alphahat_chain[i][j][k]);
      }
  cout << "\n";
  cout.flush();

  fclose(fhnd);


  gsl_matrix_free(W);


  return 0;

}
