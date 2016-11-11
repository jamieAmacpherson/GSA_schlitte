/*=============================================================================
schlitter.c: Application of the Schlitter formula to MD trajectories
Copyright (C) 2016 Jens Kleinjung, Jamie MacPherson, Franca Fraternali
=============================================================================*/

#include "config.h"
#include "schlitter.h"

/*____________________________________________________________________________*/
/* compute entropy via Schlitter formula */
/* S >= S' = 1/2 k_B ln det[1 + (k_B T e^2 / h_bar^2) M \sigma] */
/* with S:entropy; k_B:Boltzmann k; T:Temperature; e:Euler e; h:Planck h;
	M:mass matrix; \sigma:covariance matrix. */ 
__inline__ static double schlitter(Eigensys *eigensys)
{
	unsigned int i;
	double S_sch = 0.;
	const double k_B = GSL_CONST_MKSA_BOLTZMANN;
	const double T = 300;
	const double e_sq = pow(M_E, 2);
	const double h_bar_sq = pow(GSL_CONST_MKSA_PLANCKS_CONSTANT_HBAR, 2);
	const double prefr = k_B * T * e_sq / h_bar_sq;

	/* assuming eigenvalues are ordered;
		skip the 6 d.o.f. of rigid body rotation/translation */
	for (i = 0; i < eigensys->eigendim - 6; ++ i) {
		S_sch += log(1 + (prefr * gsl_vector_get(eigensys->eigenval, i)));
	} 

	S_sch *= 0.5 * k_B;
}

/*____________________________________________________________________________*/
int main(int argc, char *argv[])
{
	int i, j;
	double S_sch; /* Schlitter entropy */

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
	if (! arg.silent) fprintf(stdout, "\nInput structure\n");
    read_structure(&arg, &argpdb, &pdb);

    /*____________________________________________________________________________*/
    /** read input trajectory */
	if (arg.trajInFileName) {
		if (! arg.silent) fprintf(stdout, "\nInput trajectory\n");
		read_gromos_traj(&traj, &arg, pdb.nAllAtom);
	} else {
		fprintf(stderr, "Need trajectory file\n");
		exit(1);
	}

    /*____________________________________________________________________________*/
	/* create trajectory GSL matrix;
		in each row: x1,y1,z1,x2,y2,z2,...
		in each column: t1,t2,... */
	/* trajectory matrix */
	if (! arg.silent) fprintf(stdout, "\nCreating trajectory matrix\n");
	gsl_matrix *A = gsl_matrix_alloc(traj.nFrame,
									(3 * traj.frame[0].nAtom));

	/* assign coordinates to trajectory matrix */
	for (i = 0; i < traj.nFrame; ++ i) {
		for (j = 0; j < traj.frame[0].nAtom; ++ j) {
			gsl_matrix_set(A, i, (3 * j),  (double)traj.frame[i].trajatom[j].pos.x);
			gsl_matrix_set(A, i, ((3 * j) + 1), (double)traj.frame[i].trajatom[j].pos.y);
			gsl_matrix_set(A, i, ((3 * j) + 2), (double)traj.frame[i].trajatom[j].pos.z);
			
		}
	} 
	fprintf(stderr, "\tinitialised atom matrix: %d frames x %d coordinates\n",
				traj.nFrame, (3 * traj.frame[0].nAtom));	

    /*____________________________________________________________________________*/
	/* compute covariance matrix */
	/* assuming rigid body rotations/translations have already been removed
			from the trajectory */
	if (! arg.silent) fprintf(stdout, "\nComputing covariance matrix, may take a minute ....\n");
	gsl_matrix *C = gsl_matrix_alloc((3 * traj.frame[0].nAtom), \
									 (3 * traj.frame[0].nAtom));
	covariance(A, C);

    /*____________________________________________________________________________*/
	/* compute eigenvalues */
	if (! arg.silent) fprintf(stdout, "\nCreating eigenvalue system\n");
	eigensystem(C, &eigensys);

    /*____________________________________________________________________________*/
	/* Schlitter entropy calculation */
	if (! arg.silent) fprintf(stdout, "\nEvaluating Schlitter entropy\n");
	S_sch = schlitter(&eigensys);

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

	/* GSL matrices */
	gsl_matrix_free(A);
	gsl_matrix_free(C);

	free(eigensys.eigenvec);
	free(eigensys.eigenval);

    /*____________________________________________________________________________*/
	/* terminate */
	fprintf(stdout, "\nClean termination\n\n");

    return 0;
}
