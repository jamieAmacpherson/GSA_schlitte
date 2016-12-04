/*==============================================================================
arg.h : parse command line arguments
Copyright (C) 2007 Jens Kleinjung
Read the COPYING file for license information.
==============================================================================*/

#ifndef ARG_H
#define ARG_H

#include <assert.h>
#include <getopt.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#include "argpdb.h"
#include "config.h"
#include "error.h"

/*____________________________________________________________________________*/
/* structures */

/* variables for commmand line arguments */
typedef struct  
{
    FILE *pdbInFile;
	char *pdbInFileName;
    FILE *trajInFile;
	char *trajInFileName;
	int silent;
	int multiModel;
    FILE *schlitterOutFile;
    char *schlitterOutFileName;
} Arg;

/*____________________________________________________________________________*/
/* prototypes */
int parse_args(int argc, char **argv, Arg *arg, Argpdb *argpdb);

#endif
