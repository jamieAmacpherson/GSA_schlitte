/*=============================================================================
POPS* : Parameter OPtimised Surface of proteins and nucleic acids
Copyright (C) 2002-2016 Franca Fraternali (program author, parametrisation)
Copyright (C) 2008-2016 Jens Kleinjung (modular C code)
Copyright (C) 2002 Luigi Cavallo (parametrisation)
Copyright (C) 2002 Kuang Lin and Valerie Hindie (translation to C)

POPS is free software and comes with ABSOLUTELY NO WARRANTY.
You are welcome to redistribute it under certain conditions.
Read the COPYING file for license and distribution details.

Fraternali, F. and van Gunsteren, W.F.
An efficient mean solvation force model for use in molecular dynamics simulations
  of proteins in aqueous solution.
  Journal of Molecular Biology 256 (1996) 939-948.

Fraternali, F. and Cavallo, L.
Parameter optimized surfaces (POPS*): analysis of key interactions and 
  conformational changes in the ribosome.
  Nucleic Acids Research 30 (2002) 2950-2960.

Cavallo, L., Kleinjung, J. and Fraternali, F.
POPS: A fast algorithm for solvent accessible surface areas at atomic 
  and residue level.
  Nucleic Acids Research 31 (2003) 3364-3366.

Kleinjung, J. and Fraternali, F.
  POPSCOMP: an automated interaction analysis of biomolecular complexes.
  Nucleic Acids Research 33 (2005) W342-W346.

=============================================================================*/

#include "config.h"
#include "schlitter.h"

/*____________________________________________________________________________*/
int main(int argc, char *argv[])
{
	int i;

#ifdef DEBUG
	printf("hello DEBUG");
#endif

    /*____________________________________________________________________________*/
	/* data structures and variables */
    Arg arg; /** data structure for command line arguments */
	//Traj traj; /** data structure for trajectory */

    /*____________________________________________________________________________*/
    /** parse command line arguments */
    parse_args(argc, &(argv[0]), &arg);

    /*____________________________________________________________________________*/
    /** read input trajectory */
	//if (arg.trajInFileName) {
	//	if (! arg.silent) fprintf(stdout, "Input trajectory\n");
	//	read_gromos_traj(&traj, &arg, pdb.nAllAtom);
	//}

    /*____________________________________________________________________________*/
	/* Schlitter entropy calculation */

    /*____________________________________________________________________________*/
	/** free memory */
	/* structure */
	/* trajectory */
	//if (arg.trajInFileName) {
	//	for (i = 0; i < (traj.nFrame + 1); ++ i)
	//		free(traj.frame[i].trajatom);
	//	free(traj.frame);
	//}

    /*____________________________________________________________________________*/
	/* terminate */
	fprintf(stdout, "\nClean termination\n\n");

    return 0;
}
