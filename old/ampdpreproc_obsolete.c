/* 
 * ampdpreproc
 * ===========
 *
 * Utility program to simple preprocession of signals later fed into AMPD
 * routine for peak counting. This program doeas smoothing and filtering
 * with various methods.
 * The input should be a file with a single float value on each line
 * corresponding to the uniformly sampled values of the signal
 *
 * Usage:
 *  ampdpreproc -f [infile] -s [sampling_rate]
 *
 * Optional arguments:
 *  -o [outfile]
 *  -v --verbose 
 */

#include <stdio.h>
#include <stdlib.h>

#include "filters.h"

static int verbose;


static struct option longoptions[] =
{
    {"infile", required_argument, NULL, 'f'},
    {NULL, 0, NULL, 0}

};

int main(int argc, char **argv){

    return 0;
}
