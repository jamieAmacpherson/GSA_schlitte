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

	/* compute coordinate mean values over trajectory time */
	/* the vector has length 3N */
	//gsl_vector *meancoor = gsl_vector_alloc(A->size1);
	/* for each coordinate */
	//for (i = 0; i < A->size1; ++ i) {
		/* average over all saved time steps */
	//	gsl_vector_view meancoor = gsl_matrix_column(A, i);
	//	gsl_stats_mean(meancoor.vector.data, 1, A->size2); 
	//}

	/* compute the covariance matrix */
	/* the vectors have length 'time steps' */
	gsl_vector *ci = gsl_vector_alloc(A->size2);
	gsl_vector *cj = gsl_vector_alloc(A->size2);

	/* for each coordinate i */
	for (i = 0; i < A->size1; ++ i) {
		/* and each other coordinate j */
		for (j = i; j < A->size2; ++ j) {
			/* column vector view: reference to memory */
			gsl_vector_view c_i = gsl_matrix_column(A, i);
			gsl_vector_view c_j = gsl_matrix_column(A, j);
			covar = gsl_stats_covariance(c_i.vector.data, 1, \
										 c_j.vector.data, 1, \
										 A->size2);
			/* assign covariance to symmetric matrix elements */
			gsl_matrix_set (C, i, j, covar);
			gsl_matrix_set (C, j, i, covar);
		}
	}

	gsl_vector_free(ci);
	gsl_vector_free(cj);
}

/*____________________________________________________________________________*/
/* compute eigensystem of input matrix C */
void eigensystem(gsl_matrix *C, Eigensys *eigensys)
{
    gsl_eigen_symmv_workspace *w; 

   /* allocate eigensystem */
	eigensys->eigendim = C->size1;
	/* allocate workspace */
    w = gsl_eigen_symmv_alloc(eigensys->eigendim);
	/* allocate vector for eigenvalues */
	eigensys->eigenval = gsl_vector_calloc(eigensys->eigendim);
	/* allocate matrix for eigenvectors */
    eigensys->eigenvec = gsl_matrix_calloc(eigensys->eigendim, eigensys->eigendim);

    /* compute eigenvectors and eigenvalues */
    gsl_eigen_symmv(C, eigensys->eigenval, eigensys->eigenvec, w);

	gsl_eigen_symmv_free(w); /* free workspace */
}

