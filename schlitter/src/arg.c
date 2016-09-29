/*==============================================================================
arg.c : parse command line arguments
Copyright (C) 2007 Jens Kleinjung
Read the COPYING file for license information.
==============================================================================*/

#include "config.h"
#include "arg.h"

#ifdef MPI
#include <mpi.h>
#endif
extern int nodes;
extern int my_rank;

/*____________________________________________________________________________*/
/** print version */
static void print_version()
{
	
	fprintf(stdout, "\n"
					" __   __   __   __   *\n"
					"|__) /  \\ |__) /__` \n"
					"|    \\__/ |    .__/\n"); 

	fprintf(stdout, "\n%s %s\n", PACKAGE, VERSION);
}

/*____________________________________________________________________________*/
/** print header */
static void print_header()
{
    fprintf(stdout, "\nPOPS* : Parameter OPtimised Surface of proteins and nucleic acids\n");
}

/*____________________________________________________________________________*/
/** print license */
static void print_license()
{
    fprintf(stdout, "\nCopyright (C) 2002-2016 Franca Fraternali (program author)\n"
			"Copyright (C) 2008-2016 Jens Kleinjung (modular C code)\n"
			"Copyright (C) 2002 Kuang Lin and Valerie Hindie (translation to C)\n"
			"Copyright (C) 2002 Luigi Cavallo (parametrisation)\n"
			"POPS* is free software and comes with ABSOLUTELY NO WARRANTY.\n"
			"You are welcome to redistribute it under certain conditions.\n"
			"Read the COPYING file for distribution details.\n\n");
}

/*____________________________________________________________________________*/
/** print citation */
static void print_citation()
{
	fprintf(stdout, "\nSchlitter entropy\n");
}

/*____________________________________________________________________________*/
/** set defaults */
static void set_defaults(Arg *arg)
{
    arg->trajInFileName = "";
}

/*____________________________________________________________________________*/
/** check input */
static void check_input(Arg *arg)
{
	if (strlen(arg->trajInFileName) == 0) {
		fprintf(stderr, "No trajectory file name given\n");
		exit(1);
	}
}

/*____________________________________________________________________________*/
/** print command line arguments */
static void print_args(Arg *arg)
{
    time_t now;
    time(&now);

    fprintf(stdout, "\ndate: %s", ctime(&now));
}

/*____________________________________________________________________________*/
/** parse command line long_options */
int parse_args(int argc, char **argv, Arg *arg)
{
	int c;
	const char usage[] = "\npops [--pdb ...] [OPTIONS ...]\n\
	 INPUT OPTIONS\n\
	   --trajIn <trajectory input>\t(mode: optional , type: char  , default: void)\n\
	 OUTPUT OPTIONS\n\
	   --trajOut <POPS output>\t(mode: optional , type: char  , default: popstraj.out)\n\
	 INFO OPTIONS\n\
	   --cite\t\t\t(mode: optional , type: no_arg, default: off)\n\
	   --version\t\t\t(mode: optional , type: no_arg, default: off)\n\
	   --help\n";

    if (argc < 2) {
		print_header();
        fprintf(stderr, "%s", usage);
		print_license();
        exit(0);
    }

    set_defaults(arg);

    /** long option definition */
    static struct option long_options[] = {
        {"trajIn", required_argument, 0, 1},
        {"trajOut", required_argument, 0, 2},
        {"cite", no_argument, 0, 10},
        {"version", no_argument, 0, 11},
        {"help", no_argument, 0, 12},
        {0, 0, 0, 0}
    };

    /** assign parameters to long options */
    while ((c = getopt_long(argc, argv, "1:2:10 11 12", long_options, NULL)) != -1) {
        switch(c) {
            case 1:
                arg->trajInFileName = optarg;
                break;
            case 2:
                arg->trajOutFileName = optarg;
                break;
            case 10:
                print_citation();
                exit(0);
            case 11:
				print_version();
				print_license();
                exit(0);
            case 12:
                fprintf(stderr, "%s", usage);
				print_license();
                exit(0);
            default:
                fprintf(stderr, "%s", usage);
				print_license();
                exit(1);
        }
    }

	check_input(arg);
    print_header();
    print_args(arg);

    return 0;
}

