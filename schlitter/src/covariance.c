/*==============================================================================
covariance.c : compute covariance of GSL matrix
Copyright (C) 2016 Jens Kleinjung
Read the COPYING file for license information.
==============================================================================*/

#include "covariance.h"

/*____________________________________________________________________________*/
/* compute covariance matrix C from input matrix A */
void covariance(gsl_matrix *A, gsl_matrix *C)
{
	unsigned int i, j;
	double covar = 0.;
	gsl_vector_view ci;
	gsl_vector_view cj;

	for (i = 0; i < A->size1; ++ i) {
		for (j = i; j < A->size2; ++ j) {
			ci = gsl_matrix_column(A, i);
			cj = gsl_matrix_column(A, j);
			covar = gsl_stats_covariance(ci.vector.data, ci.vector.stride, \
										 cj.vector.data, cj.vector.stride, \
										 ci.vector.size);
			gsl_matrix_set (C, i, j, covar);
		}
	}
}

/*____________________________________________________________________________*/
/* compute eigensystem of input matrix C */
void eigensystem(gsl_matrix *C, Eigensys *eigensys)
{
    gsl_eigen_symmv_workspace *w; 

   /* allocate eigensystem */
	/* allocate workspace for 3x3 matrix */
    w = gsl_eigen_symmv_alloc(C->size1);
	/* allocate vector for eigenvalues */
	eigensys->eigenval = gsl_vector_calloc(C->size1);
	/* allocate matrix for eigenvectors */
    eigensys->eigenvec = gsl_matrix_calloc(C->size1, C->size1);
	
}

