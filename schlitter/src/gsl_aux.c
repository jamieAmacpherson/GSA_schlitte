/*=============================================================================
gsl_aux.c: Auxiliary functions for GSL 
Copyright (C) 2016 Jens Kleinjung, Jamie MacPherson, Franca Fraternali
=============================================================================*/

#include "gsl_aux.h"

/*____________________________________________________________________________*/
/* print formatted GSL vector */
int printf_gsl_vector(FILE *f, const gsl_vector *v)
{
	int status, n = 0;

   	for (size_t i = 0; i < v->size; i++) {
		for (size_t j = 0; j < 1; j++) {
			if ((status = fprintf(f, "%g ", gsl_vector_get(v, i))) < 0)
				return -1;
			n += status;
		}

		if ((status = fprintf(f, "\n")) < 0)
			return -1;
		n += status;
	}

	return n;
}

/*____________________________________________________________________________*/
/* print formatted GSL matrix */
/* taken from https://gist.github.com/jmbr/668067 */
int printf_gsl_matrix(FILE *f, const gsl_matrix *m)
{
	int status, n = 0;

   	for (size_t i = 0; i < m->size1; i++) {
		for (size_t j = 0; j < m->size2; j++) {
			if ((status = fprintf(f, "%g ", gsl_matrix_get(m, i, j))) < 0)
				return -1;
			n += status;
		}

		if ((status = fprintf(f, "\n")) < 0)
			return -1;
		n += status;
	}

	return n;
}

