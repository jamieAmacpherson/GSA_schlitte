/*==============================================================================
arg.c : parse command line arguments
Copyright (C) 2007 Jens Kleinjung
Read the COPYING file for license information.
==============================================================================*/

#include "config.h"
#include "arg.h"

/*____________________________________________________________________________*/
/** print version */
static void print_version()
{
	fprintf(stdout, "\nschlitter\n");
	fprintf(stdout, "\n%s %s\n", PACKAGE, VERSION);
}

/*____________________________________________________________________________*/
/** print header */
static void print_header()
{
    fprintf(stdout, "\nschlitter: Application of the Schlitter formula to trajectories\n");
}

/*____________________________________________________________________________*/
/** print license */
static void print_license()
{
    fprintf(stdout, "\nCopyright (C) 2016 Jens Kleinjung, Jamie MacPherson, Franca Fraternali \n"
			"'schlitter' is free software and comes with ABSOLUTELY NO WARRANTY.\n"
			"You are welcome to redistribute it under certain conditions.\n"
			"Read the COPYING file for distribution details.\n\n");
}

/*____________________________________________________________________________*/
/** print citation */
static void print_citation()
{
	fprintf(stdout, "\nnone yet\n");
}

/*____________________________________________________________________________*/
/** set defaults */
static void set_defaults(Arg *arg, Argpdb *argpdb)
{
    arg->pdbInFileName = "";
	arg->trajInFileName = 0; /* trajectory input file */ 
	argpdb->coarse = 0; /* Calpha-only computation [0,1] */
	argpdb->hydrogens = 0; /* hydrogens [0,1] */
	arg->silent = 0; /* suppress stdout */
    arg->schlitterOutFileName = "schlitter.out";
}

/*____________________________________________________________________________*/
/** check input */
static void check_input(Arg *arg, Argpdb *argpdb)
{
	if (strlen(arg->pdbInFileName) == 0)
		Error("Invalid PDB file name");
	assert(argpdb->coarse == 0 || argpdb->coarse == 1);
    assert(argpdb->hydrogens == 0 || argpdb->hydrogens == 1);
	assert(strlen(arg->schlitterOutFileName) > 0);
}

/*____________________________________________________________________________*/
/** print command line arguments */
static void print_args(Arg *arg, Argpdb *argpdb)
{
    time_t now;
    time(&now);

    fprintf(stdout, "\ndate: %s", ctime(&now));
    fprintf(stdout, \
                    "pdb: %s\n"
                    "traj: %s\n"
                    "coarse: %d\n",
        arg->pdbInFileName, arg->trajInFileName, argpdb->coarse);
    fflush(stdout);
}

/*____________________________________________________________________________*/
/** parse command line long_options */
int parse_args(int argc, char **argv, Arg *arg, Argpdb *argpdb)
{
	int c;
	const char usage[] = "\nschlitter [--pdb ...] [--traj ...] [OPTIONS ...]\n\
	 INPUT OPTIONS\n\
	   --pdb <PDB input>\t\t(mode: mandatory, type: char  , default: void)\n\
	   --traj <trajectory input>\t(mode: mandatory , type: char  , default: void)\n\
	 MODE OPTIONS\n\
	   --coarse\t\t\t(mode: optional , type: no_arg, default: off)\n\
	   --silent\t\t\t(mode: optional , type: no_arg, default: off)\n\
	 OUTPUT OPTIONS\n\
	   --schlitterOut <output>\t(mode: optional , type: char  , default: schlitter.out)\n\
	 INFO OPTIONS\n\
	   --cite\t\t\t(mode: optional , type: no_arg, default: off)\n\
	   --version\t\t\t(mode: optional , type: no_arg, default: off)\n\
	   --help\n";

    if (argc < 3) {
		print_header();
        fprintf(stderr, "%s", usage);
		print_license();
        exit(0);
    }

    set_defaults(arg, argpdb);

    /** long option definition */
    static struct option long_options[] = {
        {"pdb", required_argument, 0, 1},
        {"traj", required_argument, 0, 2},
        {"coarse", no_argument, 0, 3},
        {"schlitterOut", required_argument, 0, 6},
        {"silent", no_argument, 0, 13},
        {"cite", no_argument, 0, 25},
        {"version", no_argument, 0, 26},
        {"help", no_argument, 0, 27},
        {0, 0, 0, 0}
    };

    /** assign parameters to long options */
    while ((c = getopt_long(argc, argv, "1:2:3 5:6:13 25 26 27", long_options, NULL)) != -1) {
        switch(c) {
            case 1:
                arg->pdbInFileName = optarg;
                break;
            case 2:
                arg->trajInFileName = optarg;
                break;
            case 3:
                argpdb->coarse = 1;
                break;
            case 6:
                arg->schlitterOutFileName = optarg;
                break;
            case 13:
                arg->silent = 1;
                break;
            case 25:
                print_citation();
                exit(0);
            case 26:
				print_version();
				print_license();
                exit(0);
            case 27:
                fprintf(stderr, "%s", usage);
				print_license();
                exit(0);
            default:
                fprintf(stderr, "%s", usage);
				print_license();
                exit(1);
        }
    }

	check_input(arg, argpdb);
    print_header();
    print_args(arg, argpdb);

    return 0;
}

