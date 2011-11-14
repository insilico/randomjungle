/*
 * mvt.cpp
 *
 *  Created on: 07.04.2009
 *      Author: schwarz
 */

#include "mvt.h"


extern"C" {
// gfortran
int gsl_rng_uniform_fortran_(double *x, gsl_rng **r, int pointerwaste) {
  *x = gsl_rng_uniform(*r);
  return(1);
}

// g77
int gsl_rng_uniform_fortran__(double *x, gsl_rng **r, int pointerwaste) {
  *x = gsl_rng_uniform(*r);
  return(1);
}
}
