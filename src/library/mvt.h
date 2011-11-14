/*
 * mvt.h
 *
 *  Created on: 07.04.2009
 *      Author: schwarz
 */

#ifndef MVT_H_
#define MVT_H_

#include <gsl/gsl_rng.h>
#include <gsl/gsl_cdf.h>

extern"C" {
int gsl_rng_uniform_fortran_(double *x, gsl_rng **r, int pointerwaste); //gfortran
int gsl_rng_uniform_fortran__(double *x, gsl_rng **r, int pointerwaste); //g77
}

extern"C" {
void mvtdst_(
  int *n, int *nu, double *lower, double *upper,
  int *infin, double *corr, double *delta,
  int *maxpts, double *abseps, double *releps, double *error,
  double *value, int *inform,
  gsl_rng **r
);
}


#endif /* MVT_H_ */
