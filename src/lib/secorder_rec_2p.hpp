#include <gsl/gsl_rng.h>
#include <gsl/gsl_matrix.h>

gsl_matrix* secorder_rec_2p(int N_nodes[2], double (*p)[2],
			    double (*alpha_recip)[2], 
			    double (*alpha_conv_hom)[2], 
			    double (*cc_conv_mixed)[2][2], 
			    double (*alpha_div_hom)[2],
			    double (*cc_div_mixed)[2][2],
			    double (*alpha_chain)[2][2], gsl_rng *r);
