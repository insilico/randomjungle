//#define LOGREGDEBUG

/*
 *            Logistic.cpp (org: logistic.c)
 *
 *  fixed/reworked by Daniel F Schwarz 2007
 *  (org: Simon Bonner 2005 / sbonner@stat.sfu.ca)
 *
 */

#include <math.h>
#include <limits>
#include <iostream>
#include <gsl/gsl_math.h>
#include <gsl/gsl_sf_exp.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_permutation.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_cblas.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_errno.h>

#include "treedefs.h"
#include "Logistic.h"
#include "Exception.h"

#ifndef MPFR_MYPREC
#define MPFR_MYPREC 512
#endif

int
Logistic::logistic (gsl_vector *z, gsl_vector *b, double &val)
{
	/* Computes the logistic function with parameters b at z:
	 * 			f(z|b) = exp(b'z)/(1+exp(b'z))
	 *                               1/1+exp(-b'z)
	 * Input: 	z = vector of covariates (starting with 1 for the intercept),
	 * 			b = vector of coefficients
	 * Output: returns logistic(z)
	 */

	int status = 0;		/* GSL Error Status */
	double eta;	/* eta=b'z, val=f(z|b) */

  // Compute the linear predictor
	status = gsl_blas_ddot (z, b, &eta);
	if (status != 0) return status;

	gsl_sf_result_e10 expon;	// expon=exp(eta)
#ifdef LOGREGDEBUG
  std::cout << "eta = " << eta << std::endl;
#endif
  gsl_error_handler_t* gsl_eh = gsl_set_error_handler_off();
  status = gsl_sf_exp_e10_e(-eta, &expon);
  gsl_set_error_handler(gsl_eh);
#ifdef LOGREGDEBUG
  std::cout << "expon.val = " << expon.val << std::endl;
#endif
  val = 1 / (1 + expon.val);

	if (status != 0) return status;

	return 0;
}

int
Logistic::comp_w_z (gsl_vector * X, double *y, gsl_vector * b, gsl_vector * w,
	  gsl_vector * z)
{
	/* Computes diagonal elements of diag. matrix W
	 * and elements of vector z in WLS equation.
	 * Current implentation for case of single covariate only
	 */
	int status = 0;		/* GSL error status */
	int n = (*X).size;	/* n=number of observations */
	int i;			/* Counter */
	double p, eta;		/* p=probability of event, eta=logit(p) */
	double tmp_double;	/* Temp. value for calculations */
	gsl_vector *tmp_vec;	/* Temp. vector for calculations */
	double squaredError = 0;

	tmp_vec = gsl_vector_alloc (2);
	gsl_vector_set(tmp_vec, 0, 1.0);

#ifdef LOGREGDEBUG
  std::cout << "n = " << n << std::endl;
  std::cout << "b1 = " << gsl_vector_get(b, 0) << std::endl;
  std::cout << "b2 = " << gsl_vector_get(b, 1) << std::endl;
#endif

  for (i = 0; i < n; i++) {
    gsl_vector_set(tmp_vec, 1, gsl_vector_get(X, i));

		status = logistic(tmp_vec, b, p);
		if (status != 0) {
      gsl_vector_free (tmp_vec);
      return status;
    }

		tmp_double = p * (1 - p); //wrong differentiation?
//		tmp_double = gsl_vector_get(b, 1) * p * (1 - p);

		status = gsl_blas_ddot (tmp_vec, b, &eta);
		if (status != 0) {
      gsl_vector_free (tmp_vec);
      return status;
    }

		gsl_vector_set (w, i, p * (1 - p));
		gsl_vector_set (z, i, eta + (y[i] - p) / tmp_double);

		squaredError += (y[i] - p) * (y[i] - p);
	}
#ifdef LOGREGDEBUG
		std::cout << "squaredError == " << squaredError << std::endl;
#endif

	/* Clean-up */
	gsl_vector_free (tmp_vec);
	return 0;
}

int
Logistic::gsl_vector_sum (gsl_vector * w, double &sum)
{
	/* Sums all elements in a gsl_vector */
	gsl_vector *ones;
	int status = 0;

	/* Create vector of ones */
	ones = gsl_vector_alloc((*w).size);
	gsl_vector_set_all(ones, 1.0);

	/* Compute dot product */
	status = gsl_blas_ddot(w, ones, &sum);
	if (status != 0) {
    gsl_vector_free(ones);
    return status;
  }

	/* Clean-up */
	gsl_vector_free(ones);

	return 0;
}

int
Logistic::left_side_matrix (gsl_matrix * inf, gsl_vector * X, gsl_vector * w)
{
	/* Computes the matrix X'WX */
	/* Current implementation for case of single covariate only */
	int status = 0;		/* GSL error status */
	int n = (*X).size;	/* Number of observations */
	double tmp_double;	/* Temporary double value for calculations */
	gsl_vector *tmp_vec = gsl_vector_alloc (n);	/* Temporary vector of length n for calculations */

	/* Element 1,1: sum of w */
	status = gsl_vector_sum(w, tmp_double);
	if (status != 0) {
    gsl_vector_free(tmp_vec);
    return status;
  }
	gsl_matrix_set (inf, 0, 0, tmp_double);

	/* Element 1,2 and 2,1: dot product of w and X */
	status = gsl_blas_ddot (w, X, &tmp_double);
	if (status != 0) {
    gsl_vector_free(tmp_vec);
    return status;
  }

	gsl_matrix_set (inf, 0, 1, tmp_double);
	gsl_matrix_set (inf, 1, 0, tmp_double);

	/* Element 2,2: dot product of w with X^2 */
	gsl_vector_memcpy (tmp_vec, X);
	status = gsl_vector_mul (tmp_vec, tmp_vec);
	if (status != 0) {
    gsl_vector_free(tmp_vec);
    return status;
  }

	status = gsl_blas_ddot (w, tmp_vec, &tmp_double);
	if (status != 0) {
    gsl_vector_free(tmp_vec);
    return status;
  }

	gsl_matrix_set (inf, 1, 1, tmp_double);

	/* Clean up */
	gsl_vector_free(tmp_vec);

  return 0;
}

int
Logistic::right_side_vector(
  gsl_vector * rs,
  gsl_vector * X,
  gsl_vector * w,
	gsl_vector * z
  ) {

	/* Computes the vector X'Wz */
	/* Current implementation for case of single covariate only */
	int status = 0;		/* GSL error status */
	int n = (*X).size;	/* Number of data points */
	double tmp_double;	/* Temporary double value for calculations */
	gsl_vector *tmp_vec = gsl_vector_alloc (n);	/* Temporary vector of length n for calculations */

	/* Element 1: dot product of w and z */
	status = gsl_blas_ddot(w, z, &tmp_double);
	if (status != 0) {
    gsl_vector_free(tmp_vec);
    return status;
  }

	gsl_vector_set(rs, 0, tmp_double);

	/* Element 2: dot product of w and X*z */
	gsl_vector_memcpy (tmp_vec, X);
	gsl_vector_mul (tmp_vec, z);
	status = gsl_blas_ddot (w, tmp_vec, &tmp_double);
	if (status != 0) {
    gsl_vector_free(tmp_vec);
    return status;
  }

	gsl_vector_set (rs, 1, tmp_double);

	/* Clean up */
	gsl_vector_free (tmp_vec);

	return 0;
}

int
Logistic::solve_sym (gsl_matrix * A, gsl_vector * x, gsl_vector * b, gsl_matrix * A_inv)
{
/* Solves the system Ax=b where A is a symmetric matrix*/
/* Returns inverse of A as by product of calculation */
/* Currently implemented only for the case of dimension 2 */
	int status = 0;		/* GSL error status */
	double a11, a12, a22;	/* Elements of symmetric matrix A */

	/* 1. Unpack A for simplicity */
	a11 = gsl_matrix_get (A, 0, 0);
	a12 = gsl_matrix_get (A, 0, 1);
	a22 = gsl_matrix_get (A, 1, 1);

	/* 2. Compute the inverse of A directly */
	/* Re-arrange elements of A */
	gsl_matrix_set (A_inv, 0, 0, a22);
	gsl_matrix_set (A_inv, 1, 0, -1 * a12);
	gsl_matrix_set (A_inv, 0, 1, -1 * a12);
	gsl_matrix_set (A_inv, 1, 1, a11);

	/* Divide by determinant of A */
	gsl_matrix_scale (A_inv, 1 / (a11 * a22 - a12 * a12));

	/* 3. Multiply the right side by the inverse information */
  status = gsl_blas_dgemv (CblasNoTrans, 1, A_inv, b, 0, x);
	if (status != 0) return status;

  return 0;
}

int
Logistic::logistic_reg_iter(
  gsl_vector * X,
  double *resp,
  gsl_vector * b,
  gsl_matrix * cov,
  double *diff
  ){

/*
 * Computes one iteration of the Fisher scoring algorithm for case of one scalar
 * covariate
 * Inputs: X=vector of predictors, resp=vector of {0,1} responses
 * Returns: b=update coefficients, cov=variance matrix, diff=max change in coef.
 */
	/*
	 * b_new = b_old + (X'WX)^-1X'W\widetilde{Y}=(X'WX)^-1X'WZ
	 * =>
	 * The new estimates, b, are computed using the weighted LS equation:
	 * X'WXb_new=X'WZ
	 *
	 * W are the weights (diag. matrix):
	 * 		$\displaystyle {\boldsymbol{W}}
	 * 		= {\text{diag}}\left( 	\frac{G'(\eta_1)^2}{V(\mu_1)},
	 * 								\ldots,
	 * 								\frac{G'(\eta_n)^2}{V(\mu_n)}\right)$
	 * Z are the adjusted dependent variables:
	 * 		Z_i = \boldsymbol{X}_i^\top \boldsymbol{\beta}^{\text{old}}
	 * 			+ \frac{Y_i-\mu_i}{G'(\eta_i)}\,.$
	 * 		with
	 * 		$ \mu_i=G(\eta_i)
	 * 		= G(\boldsymbol{X}_i^\top \boldsymbol{\beta}) =b'(\theta_i)$
	 * X are the explanatory variables
	 *
	 *
	 *
	 */

	int status = 0;		/* GSL error status */
	int n = (*X).size;	/* n=number of observations */
	gsl_vector *w = gsl_vector_alloc (n);	/* w=p*(1-p) for each observation */
	gsl_matrix *ls = gsl_matrix_alloc (2, 2);	/* Matrix on left side of WLS equation: X'WX */
	gsl_vector *rs = gsl_vector_alloc (2);	/* Vector on right side of WLS equation: X'Wz */
	gsl_vector *z = gsl_vector_calloc (n);	/* Vector z in WLS eqtuation */
	gsl_vector *b_old = gsl_vector_alloc (2);	/* Stores coeff. from previous iteration */

/* Store current coefficient values in b_old */
	status = gsl_vector_memcpy (b_old, b);
	if (status != 0) {
    gsl_vector_free (w);
    gsl_vector_free (rs);
    gsl_matrix_free (ls);
    gsl_vector_free (z);
    gsl_vector_free (b_old);
    return status;
  }


/* Compute values of w and z for each data point  */
	status = comp_w_z (X, resp, b_old, w, z);
	if (status != 0) {
    gsl_vector_free (w);
    gsl_vector_free (rs);
    gsl_matrix_free (ls);
    gsl_vector_free (z);
    gsl_vector_free (b_old);
    return status;
  }

/* Compute the matrix on the left side of the equation: X'WX */
	status = left_side_matrix (ls, X, w);
	if (status != 0) {
    gsl_vector_free (w);
    gsl_vector_free (rs);
    gsl_matrix_free (ls);
    gsl_vector_free (z);
    gsl_vector_free (b_old);
    return status;
  }

/* Compute vector on right side of equation: X'Wz */
	status = right_side_vector (rs, X, w, z);
	if (status != 0) {
    gsl_vector_free (w);
    gsl_vector_free (rs);
    gsl_matrix_free (ls);
    gsl_vector_free (z);
    gsl_vector_free (b_old);
    return status;
  }

/* Solve the Fisher scoring equation: X'WX b=X'Wz */
	/* The updated variance-covariance matrix is a by product of this calculation */
	status = solve_sym (ls, b, rs, cov);
	if (status != 0) {
    gsl_vector_free (w);
    gsl_vector_free (rs);
    gsl_matrix_free (ls);
    gsl_vector_free (z);
    gsl_vector_free (b_old);
    return status;
  }

/* Compute the maximum difference between the old and new coefficients to test for convergence */
	status = gsl_vector_sub (b_old, b);
	if (status != 0) {
    gsl_vector_free (w);
    gsl_vector_free (rs);
    gsl_matrix_free (ls);
    gsl_vector_free (z);
    gsl_vector_free (b_old);
    return status;
  }

	*diff = GSL_MAX (fabs (gsl_vector_get (b_old, 0)),
			 fabs (gsl_vector_get (b_old, 1)));

/* Clean-up */

	gsl_vector_free (w);
	gsl_vector_free (rs);
	gsl_matrix_free (ls);
	gsl_vector_free (z);
	gsl_vector_free (b_old);

	return 0;
}


/*
 * Conducts a logistic regression of resp (y) on predictor (X)
 * Estimates are computed using the Fisher scoring algorithm
 * Starting values are the vector b provided
 * Inputs: X=design matrix, resp=response vector, epsilon=tolerance for converg.
 * Returns: b=estimated parameters, cov=variance-covariance matrix
 */
void
Logistic::logistic_reg(//INPUT
						gsl_vector* X,
						double* resp,
						//OUTPUT
						gsl_vector* b,
						gsl_matrix* cov,
						//INPUT
						double epsilon,
						bool &isNotConverging) {

	// Keeps track of the maximum difference on each iteration
	double diff, lastDiff;
	uli_t steps = 0;
	uli_t stepStop = 25;
	double meanDiff = 0;
	double cutoffConv = 0;
	gsl_vector* bold = gsl_vector_alloc(2);
	int status;

	/* While diff > epsilon iterate fisher scoring algorithm */
	diff = 2 * epsilon;

	while (steps < stepStop) {
    status = gsl_vector_memcpy(bold, b);
    if (status != 0) throw Exception(gsl_strerror(status));

    // calc cutoff as a convegency criteria
		cutoffConv = -gsl_vector_get(b, 0) / gsl_vector_get(b, 1);

		// iterate one step in logistic regression (via fisher information)
		status = logistic_reg_iter (X, resp, b, cov, &diff);
		if (status != 0) {
		  isNotConverging = true;
		  return;
      //throw Exception(gsl_strerror(status));
    }

		lastDiff = meanDiff;
		meanDiff += diff;

    if (isnan(gsl_vector_get(b, 0)) || isnan(gsl_vector_get(b, 0))) {
      gsl_vector_memcpy(b, bold);
      break;
    }

#ifdef LOGREGDEBUG
		std::cout 	<< "b0 " << gsl_vector_get(b, 0) << std::endl;
		std::cout 	<< "b1 " << gsl_vector_get(b, 1) << std::endl;
		std::cout 	<< "meanDiff/(1 + steps) == "
					<< meanDiff/(1 + steps) << std::endl;
		std::cout 	<< "lastDiff/(steps) == "
					<< lastDiff/steps << std::endl;
		std::cout 	<< "diff == "
					<< diff << std::endl;
		std::cout 	<< "cutoffConv == "
					<< cutoffConv << std::endl;
#endif

		++steps;
	}
#ifdef LOGREGDEBUG
	std::cout 	<< "error == " <<
	fabs(-gsl_vector_get(b, 0) / gsl_vector_get(b, 1)
				- cutoffConv)
				<< std::endl;
#endif
	if (fabs(-gsl_vector_get(b, 0) / gsl_vector_get(b, 1)
			- cutoffConv) > epsilon) {
		isNotConverging = true;
	} else {
		isNotConverging = false;
	}

	gsl_vector_free(bold);

/*
	if (fabs(meanDiff/steps-lastDiff/(steps-1)) > epsilon ) {
#ifdef LOGREGDEBUG
		std::cout 	<< "LR is not converging." << std::endl;
#endif
		isNotConverging = true;
	} else {
		isNotConverging = false;
	}
#ifdef LOGREGDEBUG
	std::cout << std::endl;
#endif
*/
}
