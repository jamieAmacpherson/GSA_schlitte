/*==============================================================================
covariance.c : compute covariance of GSL matrix
Copyright (C) 2016 Jens Kleinjung
Read the COPYING file for license information.
==============================================================================*/

#include "covariance.h"

/*____________________________________________________________________________*/
/* compute covariance matrix C from input matrix A */
/* using the following R implementation as example:
	1/(NROW(A)-1) * crossprod(scale(A, TRUE, FALSE)) */
/* 1. centre the matrix 'A' columns by subtracting the column mean;
   2. compute the cross product 'A x A';
   3. normalise with the 1/(n-1) prefactor. */
void cov(gsl_matrix *A, gsl_matrix *C)
{
	unsigned int i, j;
	double colsum, colmean;
	/* 1. centre the matrix 'A' columns by subtracting the column mean */
	/* 'Ac': centred matrix 'A' */
	gsl_matrix *Ac = gsl_matrix_alloc(A->size1, A->size2);
	/* column vector has length 'time steps' */
	gsl_vector *cj = gsl_vector_alloc(A->size2);
	/* for all columns of 'A' */
	for (j = 0; j < A->size2; ++ j) {
		gsl_vector_view c_j = gsl_matrix_column(A, j);
		/* compute column mean */
		for (i = 0, colsum = 0; i < A->size1; ++ i) {
			colsum += gsl_vector_get(&(c_j.vector), i);
		}
		colmean = colsum / A->size1;

		/* assign centred column A[ ,j] to Ac[ ,j] */
		for (i = 0; i < A->size1; ++ i) {
			gsl_matrix_set(Ac, i, j, gsl_matrix_get(A, i, j) - colmean);
		}
	}
#ifdef DEBUG
	FILE *ceOut = safe_open("tmp_centre_C.dat", "w");
	printf_gsl_matrix(ceOut, Ac);
	fclose(ceOut);
#endif

	/* 2. compute the covariance matrix 'C' from the cross-product:
		C = A x A */
	/* the crossproduct 'DGEMM' in BLAS:
		C := alpha * A x B + beta * C */
	gsl_blas_dgemm(CblasConjTrans, CblasNoTrans, 1., Ac, Ac, 0., C);

#ifdef DEBUG
	FILE *croOut = safe_open("tmp_crossprod_C.dat", "w");
	printf_gsl_matrix(croOut, C);
	fclose(croOut);
#endif

	/* 3. normalise cocariance matrix 'C' with the 1/(NROW(A)-1) prefactor. */
	const double x = (double)1. / (Ac->size1 - (double)1.);
	gsl_matrix_scale(C, x);

#ifdef DEBUG
	fprintf(stderr, "\tcovariance scaling: %lf\n", x);
	FILE *coOut = safe_open("tmp_covar_C.dat", "w");
	printf_gsl_matrix(coOut, C);
	fclose(coOut);
#endif
}

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

