/*
 * ampdr.c
 *
 * Implementation of AMPD for CPU with some optimization .
 *
 */

#include <math.h>
#include "ampdr.h"

/* Constants are picked according to imput data type.
 *
 * Respiration=0
 * Pulsoxymetry=1
 * ECG=2
 *
 */

#define DATATPYE 1


/**
 * Main routine for peak detection on a dataseries
 *
 * @param data      Input dataseries, previosly loaded from file
 * @param n         Length of dataseries
 * @param lms       Local maxima scalogram matrix
 * @param rlm       Rscaled local maxima scalogram
 * @param gamma     Vector coming from the row-wise summation of the LMS matrix.
 * @param sigma     Vector from the column-wise standard deviation of the rescaled
 *                  LMS matrix
 * @param peaks     Vector containing the indices of peaks corresponding to the 
 *                  original input dataseries
 *
 * @return          Number of peaks if successful, -1 on error.
 */

int ampdcpu(float *data,int n,struct fmtx *lms,float *gam,float *sig, int *pks){

    int n_pks;

    return n_pks;


}
