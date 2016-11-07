/*=============================================================================
schlitter.c: Application of the Schlitter formula to MD trajectories
Copyright (C) 2016 Jens Kleinjung, Jamie MacPherson, Franca Fraternali
=============================================================================*/

#include "config.h"
#include "schlitter.h"

/*____________________________________________________________________________*/
int main(int argc, char *argv[])
{
	int i;

	gsl_matrix *A = gsl_matrix_alloc (100, 100);
	gsl_matrix *C = gsl_matrix_alloc (100, 100);

    /*____________________________________________________________________________*/
	/* data structures and variables */
    Arg arg; /** data structure for command line arguments */
    Argpdb argpdb; /** data structure for PDB command line arguments */
	Str pdb; /** data structure for input PDB molecule */
	Traj traj; /** data structure for trajectory */
	Eigensys eigensys; /* eigen system */

    /*____________________________________________________________________________*/
    /** parse command line arguments */
    parse_args(argc, &(argv[0]), &arg, &argpdb);

    /*____________________________________________________________________________*/
    /** read input structure */
	if (! arg.silent) fprintf(stdout, "Input structure\n");
    read_structure(&arg, &argpdb, &pdb);

    /*____________________________________________________________________________*/
    /** read input trajectory */
	if (arg.trajInFileName) {
		if (! arg.silent) fprintf(stdout, "Input trajectory\n");
		read_gromos_traj(&traj, &arg, pdb.nAllAtom);
	}

    /*____________________________________________________________________________*/
	/* compute covariance matrix */
    for (i = 0; i < traj.nFrame; ++ i) {
		if (! arg.silent) {
            (((i+1) % 50) != 0) ? fprintf(stdout, ".") : fprintf(stdout, "%d\n\t", (i + 1));
			fflush(stdout);
        }
        assert(traj.frame[i].nAtom == pdb.nAllAtom);
	}

    /*____________________________________________________________________________*/
	/* compute eigenvalues */
	eigensystem(C, &eigensys);

    /*____________________________________________________________________________*/
	/* Schlitter entropy calculation */
	/* S >= S' = 1/2 k_B ln det[1 + (k_B T e^2 / h_bar^2) M \sigma] */
	/* with S:entropy; k_B:Boltzmann k; T:Temperature; e:Euler e; h:Planck h;
		M:mass matrix; \sigma:covariance matrix. */ 


    /*____________________________________________________________________________*/
	/** free memory */
	/* structure */
    free(pdb.sequence.name);
    free(pdb.atom);
    free(pdb.atomMap);
    free(pdb.sequence.res);
	/* trajectory */
	if (arg.trajInFileName) {
		for (i = 0; i < (traj.nFrame + 1); ++ i)
			free(traj.frame[i].trajatom);
		free(traj.frame);
	}

    /*____________________________________________________________________________*/
	/* terminate */
	fprintf(stdout, "\nClean termination\n\n");

    return 0;
}
