/*=============================================================================
gsl_aux.c: Auxiliary functions for GSL 

Copyright (C) 2016 Jens Kleinjung, Jamie MacPherson, Franca Fraternali
=============================================================================*/

#ifndef GSL_AUX_H
#define GSL_AUX_H

#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <stdlib.h>
#include <stdio.h>
#include "config.h"

int printf_gsl_vector(FILE *f, const gsl_vector *v);
int printf_gsl_matrix(FILE *f, const gsl_matrix *m);

#endif

