/*==============================================================================
covariance.h : compute covariance of GSL matrix
Copyright (C) 2016 Jens Kleinjung
Read the COPYING file for license information.
==============================================================================*/

#ifndef COVARIANCE_H
#define COVARIANCE_H

#include <gsl/gsl_blas.h>
#include <gsl/gsl_eigen.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_statistics_double.h>
#include "config.h"
#include "gsl_aux.h"
#include "safe.h"

typedef struct
{
	unsigned int eigendim;
    gsl_vector *eigenval;
    gsl_matrix *eigenvec;
} Eigensys;

void cov(gsl_matrix *A, gsl_matrix *C);
void covariance(gsl_matrix *A, gsl_matrix *C);
void eigensystem(gsl_matrix *C, Eigensys *eigensys);

#endif

