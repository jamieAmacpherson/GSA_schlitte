/*===============================================================================
getgromacs.h : Read GROMOS96 trajectory file
Copyright (C) 2016 Jens Kleinjung
Read the COPYING file for license information.
================================================================================*/

#ifndef GETGROMACS_H
#define GETGROMACS_H

#include <math.h>
#include <regex.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "arg.h"
#include "config.h"
#include "safe.h"

int read_gromos_traj(Arg *arg, int protEnd);

#endif
