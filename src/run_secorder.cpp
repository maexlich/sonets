#include <iostream>
#include <string.h>
#include <fcntl.h>
#include <sys/stat.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_math.h>
#include "lib/secorder_rec_1p.hpp"
#include "lib/calc_stats_1p.hpp"
#include "lib/calc_rhos.hpp"


using namespace std;

int main(int argc, char *argv[]) {
  

  ///////////////////////////////////////////////////////////////
  // Read seven input parameters 
  // N_nodes      number of nodes in the network
  // p            the probability of any one given connection
  // alpha_recip  determines covariance of reciprocal connections
  // alpha_conv   determines covariance of convergent connections
  // alpha_div    determines covariance of divergent connections
  // cc_chain     determines covariance of chain, 
  //              relative to alpha_conv and alpha_div
  // seed         seed for the random number generator (optional)
  ///////////////////////////////////////////////////////////////

 
  if(argc < 7 || argc > 8) {
    cerr << "Requires six or seven parameters: N_nodes p alpha_recip alpha_conv alpha_div cc_chain [seed]\n";
    exit(-1);
  }

  int N_nodes = atoi(argv[1]);
  double p = atof(argv[2]);
  double alpha_recip = atof(argv[3]);
  double alpha_conv = atof(argv[4]);
  double alpha_div = atof(argv[5]);
  double cc_chain = atof(argv[6]);

  gsl_rng *rng = gsl_rng_alloc(gsl_rng_mt19937);


  int deterministic_seed=0;
  int rng_seed = 0;
  if(argc >7) {
    deterministic_seed=1;
    rng_seed = atoi(argv[7]);    
  }

  // set the seed
  // if deterministic_seed, use rng_seed for seed
  if(deterministic_seed)
    gsl_rng_set(rng,rng_seed); 
  // else use time in seconds for the seed
  else
    gsl_rng_set(rng, time(NULL));
  
  gsl_matrix *W;
  W=secorder_rec_1p(N_nodes,p,  alpha_recip, alpha_conv, alpha_div, cc_chain,
		    rng);
  
    // if failed to generate a matrix, write error and quit
  if(!W) {
    cerr << "Failed to generate the matrix\n";
    return -1;
  }

  // matrix files will be stored in data directory
  // with filename determined by the six or seven input parameters
  mkdir("data",0755);  // make data directory, if it doesn't exist
  char FNbase[200];
  if(deterministic_seed)
    sprintf(FNbase,"_%i_%1.3f_%1.3f_%1.3f_%1.3f_%1.3f_%i",N_nodes,p,alpha_recip, alpha_conv, alpha_div, cc_chain, rng_seed);
  else
    sprintf(FNbase,"_%i_%1.3f_%1.3f_%1.3f_%1.3f_%1.3f",N_nodes,p,alpha_recip, alpha_conv, alpha_div, cc_chain);
  
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

  for(int i=0; i<N_nodes;i++) {
    for(int j=0; j<N_nodes; j++) {
      fprintf(fhnd, "%i ", (int) gsl_matrix_get(W,i,j));
    }
    fprintf(fhnd,"\n");
  }
  fclose(fhnd);

  

  ////////////////////////////////////////////////////////////
  // Calculate the covariance structure the adjacency matrix
  // This should approximately agree with the input alphas
  ////////////////////////////////////////////////////////////

  cout << "Testing statistics of W...";
  cout.flush();
  double phat, alphahat_recip, alphahat_conv, alphahat_div,
    alphahat_chain, alphahat_other;
  
  calc_phat_alphahat_1p(W, N_nodes, phat, alphahat_recip, 
			alphahat_conv, alphahat_div,
			alphahat_chain, alphahat_other);

  // calculate alpha_chain from alpha_conv, alpha_div, and cc_chain
  int status;
  double rho_conv = calc_rho_given_alpha(p, p, alpha_conv, status);
  double rho_div = calc_rho_given_alpha(p, p, alpha_div, status);
  double rho_chain = cc_chain*sqrt(rho_conv*rho_div);
  double alpha_chain = calc_alpha_given_rho(p, p, rho_chain, status);

  strcpy(FN, "data/stats");
  strcat(FN, FNbase);
  strcat(FN, ".dat");
  fhnd = fopen(FN, "w");
  if(fhnd==NULL) {
    cerr << "Couldn't open outfile file " << FN << "\n";
    exit(-1);
  }
  fprintf(fhnd, "%e %e %e %e %e %e\n", phat, alphahat_recip, 
	  alphahat_conv, alphahat_div, alphahat_chain, alphahat_other);
  fprintf(fhnd, "%e %e %e %e %e %e\n", p, alpha_recip, 
	  alpha_conv, alpha_div, alpha_chain, cc_chain);
  fclose(fhnd);


  cout << "done\n";
  cout << "\nComparison of prescribed to actual statistics of matrix:\n";
  cout << "p = " << p << ", phat = " << phat << "\n";
  cout << "alphas:\n";
  cout << "alpha_recip = " << alpha_recip
       << ", alpha_conv = " << alpha_conv
       << ", alpha_div = " << alpha_div
       << ", alpha_chain = " << alpha_chain
       << ", alpha_other = 0\n";
  cout << "alphahats:\n";
  cout << "alpha_recip = " << alphahat_recip
       << ", alpha_conv = " << alphahat_conv
       << ", alpha_div = " << alphahat_div
       << ", alpha_chain = " << alphahat_chain
       << ", alpha_other = " << alphahat_other << "\n";
  
}
