#ifndef LOGISTIC_H_
#define LOGISTIC_H_

#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>



class Logistic
{
public:
	static int logistic (gsl_vector *z, gsl_vector * b, double &val);

	static void logisticreg_iter_2(gsl_matrix *X,double *resp,gsl_vector *b,gsl_matrix *inf,double *diff);

	static int comp_w_z(gsl_vector *X, double *y, gsl_vector *b, gsl_vector *w, gsl_vector *z);

	static int gsl_vector_sum(gsl_vector *w, double &sum);

	static int left_side_matrix(gsl_matrix *inf,gsl_vector *X, gsl_vector *w);

	static int right_side_vector(gsl_vector *rs,gsl_vector *X, gsl_vector *w, gsl_vector *z);

	static int solve_sym(gsl_matrix *A,gsl_vector *b, gsl_vector *x, gsl_matrix *A_inv);

	static int logistic_reg_iter (gsl_vector *X, double *resp, gsl_vector * b,gsl_matrix * cov, double *diff);

	static void logistic_reg(gsl_vector *X,double *resp, gsl_vector *b,gsl_matrix *cov, double epsilon, bool &isNotConverging);
};

#endif /*LOGISTIC_H_*/


