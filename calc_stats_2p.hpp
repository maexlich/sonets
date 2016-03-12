#include <gsl/gsl_matrix.h>

int calc_N_sec(gsl_matrix *W, int N_nodes[2], 
	       double (*N_edge)[2],
	       double (*N_recip)[2], 
	       double (*N_conv)[2][2], double (*N_div)[2][2],
	       double (*N_chain)[2][2]);
int calc_phat_alphahat(gsl_matrix *W, int N_nodes[2], 
		       double (*phat)[2],
		       double (*alphahat_recip)[2], 
		       double (*alphahat_conv)[2][2], 
		       double (*alphahat_div)[2][2],
		       double (*alphahat_chain)[2][2]);
int calc_gaus_covs(gsl_matrix *W, int N_nodes[2], 
		   double (*sigma)[2],
		   double (*cov_recip)[2], 
		   double (*cov_conv)[2][2], double (*cov_div)[2][2],
		   double (*cov_chain)[2][2],
		   double &cov_other);
