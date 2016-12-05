/*=============================================================================
schlitter.c: Application of the Schlitter formula to MD trajectories
Copyright (C) 2016 Jens Kleinjung, Jamie MacPherson, Franca Fraternali
=============================================================================*/

#include "config.h"
#include "schlitter.h"

/*____________________________________________________________________________*/
/* compute entropy via Schlitter formula */
/* S >= S' = 1/2 k_B ln det[1 + (k_B T e^2 / h_bar^2) M covar] */
/* with S:entropy; k_B:Boltzmann k; T:Temperature; e:Euler e; h:Planck h;
	M:mass vector; covar:covariance matrix, diagonalised. */
/* units: the prefactor ('prefr') units cancel out with the mass vector (kg)
	and the diagonalised distance variance matrix (m^2). Therefore the term
	in the 'log' expression is unitless and the k_B prefactor assigns the
	entropy unit J K^{-1}. */
/* more units: The atomic mass is that of a unified C$^{\alpha}$H:
	13.01864 au converted to kg.
	The distance is given is nm units and converted to m. */
__inline__ static double schlitter(Eigensys *eigensys, int nFrame)
{
	unsigned int i;
	FILE *outFile;
	FILE *csoutFile;
	double S_sch = 0.; /* Schlitter entropy per eigenvalue */
	double S_sch_cumsum = 0.; /* cumulative sum of Schlitter entropy */
	const double k_B = (double)GSL_CONST_MKSA_BOLTZMANN; /* J K^{-1} */
	const double T = 300; /* K */
	const double e_sq = pow(M_E, 2);
	const double h_bar_sq = pow((double)GSL_CONST_MKSA_PLANCKS_CONSTANT_HBAR, 2); /* J s */
	const double prefr = k_B * T * e_sq / h_bar_sq; /* kg^{-1} m^{-2} */
	const double m_CH = 13.01864; /* mass of unified C$^{\alpha}$H in au units */
	/* unified atomic mass: conversion factor from au to kg */
	const double cf_au_kg = (double)GSL_CONST_MKSA_UNIFIED_ATOMIC_MASS;
	const double cf_nmsq_msq = 1e-18; /* conversion factor from nm^2 to m^2 */
	const double Nav = (double)GSL_CONST_NUM_AVOGADRO;
	double ev = 0.;
	double logterm;

	outFile = safe_open("S_sch_C.dat", "w");
	csoutFile = safe_open("S_sch_cumsum_C.dat", "w");

	/* assuming eigenvalues are ordered;
		Skip the 6 d.o.f. of rigid body rotation/translation
		plus any whose rank is >= the number of time points.
		For example, 100 time points can only paramtrise 100 coordinates */
	
	logterm = prefr * cf_au_kg * m_CH * cf_nmsq_msq;
	for (i = 0, S_sch = 0, S_sch_cumsum = 0; (i < eigensys->eigendim - 6) && (i < nFrame - 6); ++ i) {
		ev = gsl_vector_get(eigensys->eigenval, i);
		if (ev > 0) {
			S_sch = Nav * 0.5 * k_B * log(1 + (logterm * ev));
#ifdef DEBUG
			fprintf(outFile, "%d\t%lf\t%lf\t%lf\t%e\n",
						i,
						logterm,
						logterm * ev,
						log(1 + (logterm * ev)),
						S_sch);
#else
			fprintf(outFile, "%d\t%e\n",
#endif

			S_sch_cumsum += S_sch;
			fprintf(csoutFile, "%d\t%e\n", i, S_sch_cumsum);
		} else {
			break;
		}
	}

	fclose(outFile);
	fclose(csoutFile);

	fprintf(stdout, "\tEntropy computation for %d time points\n",
				GSL_MIN((eigensys->eigendim - 6), (nFrame - 6)));

	return S_sch_cumsum;
}

/*____________________________________________________________________________*/
int main(int argc, char *argv[])
{
	int i, j;
	double S_sch; /* Schlitter entropy */
	FILE *trajOut; /* output trajectory matrix */
	FILE *covarOut; /* output covariance matrix */
	FILE *eigenvalOut; /* output eigen values */
	FILE *eigenvecOut; /* output eigen vectors */

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
	/* create trajectory GSL matrix */
	/* matrix order: time along row dimension, coordinates along column direction: */
	/* x1(t1) y1(t1) z1(t1) x2(t1) y2(t1) z2(t1) ... zn(t1)
	   x1(t2) y1(t2) z1(t2) x2(t2) y2(t2) z2(t2) ... zn(t2)
	   ...
	   x1(tn) y1(tn) z1(tn) x2(tn) y2(tn) z2(tn) ... zn(tn) */
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
	fprintf(stdout, "\tinitialised atom matrix: %d frames x %d coordinates\n",
				traj.nFrame, (3 * traj.frame[0].nAtom));	

#ifdef DEBUG
	trajOut = safe_open("debug_traj_C.dat", "w");
	printf_gsl_matrix(trajOut, A);
	fclose(trajOut);
#endif

    /*____________________________________________________________________________*/
	/* compute covariance matrix */
	/* assuming rigid body rotations/translations have already been removed
			from the trajectory */
	if (! arg.silent) fprintf(stdout, "\nComputing covariance matrix, may take a minute ....\n");
	gsl_matrix *C = gsl_matrix_alloc((3 * traj.frame[0].nAtom), \
									 (3 * traj.frame[0].nAtom));
	/* matrix C holds the covariance of matrix A */
	/* matrix order of A: see above */
	/* matrix order of C: square with coordinates in row and matrix dimension */
	cov(A, C);

    /*____________________________________________________________________________*/
	/* compute eigenvalues */
	if (! arg.silent) fprintf(stdout, "\nCreating eigenvalue system\n");

	eigensystem(C, &eigensys);

#ifdef DEBUG
	eigenvalOut = safe_open("debug_eigenval_C.dat", "w");
	printf_gsl_vector(eigenvalOut, eigensys.eigenval);
    //gsl_vector_fprintf (eigenvalOut, eigensys.eigenval, "%f");
	fclose(eigenvalOut);

	eigenvecOut = safe_open("debug_eigenvec_C.dat", "w");
	printf_gsl_matrix(eigenvecOut, eigensys.eigenvec);
    //gsl_matrix_fprintf(eigenvecOut, eigensys.eigenvec, "%f");
	fclose(eigenvecOut);
#endif

    /*____________________________________________________________________________*/
	/* Schlitter entropy calculation */
	if (! arg.silent) fprintf(stdout, "\nEvaluating Schlitter entropy\n");
	S_sch = schlitter(&eigensys, traj.nFrame);

	fprintf(stdout, "Schlitter entropy: %e J / (mol K)\n", S_sch);

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
