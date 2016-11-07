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

#ifdef DEBUG
	printf("hello DEBUG");
#endif

    /*____________________________________________________________________________*/
	/* data structures and variables */
    Arg arg; /** data structure for command line arguments */
    Argpdb argpdb; /** data structure for PDB command line arguments */
	Str pdb; /** data structure for input PDB molecule */
	Traj traj; /** data structure for trajectory */

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
