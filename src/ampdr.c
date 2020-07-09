/*
 * ampdr.c
 *
 * Implementation of AMPD for CPU with some optimization .
 *
 * Filter should be another individual function
 */

#include <stdio.h>
#include <math.h>
#include <string.h>
#include "ampdr.h"

/* Constants for different data sources. Best chioces depend on signal strength
 * data collection rate, average noise, etc.
 *
 * RESP - Respiration dataset
 * PUSL - Pulsoxymeter dataset
 * ECG - Electrocardiogram
 * DEF - default, or not specified
 */

#define RESP_TOL 0.13
#define RESP_TIMERES 0.0002
#define RESP_TIMETHRESH 0.1
//TODO
#define PULS_TOL 0.13
#define PULS_TIMERES 0.0002
#define PULS_TIMETHRESH 0.1

#define ECG_TOL 0.13
#define ECG_TIMERES 0.0002
#define ECG_TIMETHRESH 0.1

#define DEF_TOL 0.13
#define DEF_TIMERES 0.0002
#define DEF_TIMETHRESH 0.1

/**
 * Main routine for peak detection on a dataseries
 *
 * @param data      Input dataseries, previosly loaded from file
 * @param n         Length of dataseries
 *
 * The following are optional inputs. These can be given so they can be saved
 * to file once ampd is done. If nullpointer is given, ampd malloc memory for
 * subsequent calculations and frees it up once completed.
 *
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

int ampdcpu(float *data, int n, struct fmtx *lms, 
            float *gam,float *sig, int *pks, char *dtype){

    /* To keep track of nullpointer inputs. 1 means nullponter input. in order:
     * lms, gam, sig, pks
     */
    int null_inputs[4] = {0,0,0,0}; 
    if(lms == NULL)
        null_inputs[0] = 1;
    if(gam == NULL)
        null_inputs[1] = 1;
    if(sig == NULL)
        null_inputs[1] = 1;
    if(pks == NULL)
        null_inputs[1] = 1;
    int n_pks;
    int i, j;
    double tol = 0.0; int min_dst = 0;
    set_constants(dtype, &tol, &min_dst);


    return n_pks;


}
/*
 * Set smoothing window halfwidth, sigma tolerance and minimum peak distance
 * based on data type.
 */
void set_constants(char *type, double *tol, int *min_dst){

    double timeres, timethresh, smoothwin;
    if(strcmp(type, "resp")==0){
        *tol = (double)RESP_TOL;
        timeres = RESP_TIMERES;
        timethresh = RESP_TIMETHRESH;
    }
    else if(strcmp(type, "puls")==0){
        *tol = (double)PULS_TOL;
        timeres = PULS_TIMERES;
        timethresh = PULS_TIMETHRESH;
    }
    else if(strcmp(type, "ecg")==0){
        *tol = (double)ECG_TOL;
        timeres = ECG_TIMERES;
        timethresh = ECG_TIMETHRESH;
    }
    // default, but not optimized
    else{
        *tol = DEF_TOL;
        timeres = DEF_TIMERES;
        timethresh = DEF_TIMETHRESH;
    }
    (*min_dst) = (int)(timethresh / timeres);
}

