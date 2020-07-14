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


//TODO move these into main

/* Constants for different data sources. Best chioces depend on signal strength
 * data collection rate, average noise, etc.
 *
 * RESP - Respiration dataset
 * PUSL - Pulsoxymeter dataset
 * ECG - Electrocardiogram
 * DEF - default, or not specified
 *
 * SIGMA_THRESHOLD - level below sigma counts as zero
 * PEAK_THRESHOLD - minimum peak distance in seconds
 */

#define RESP_SAMPLING_RATE 200
#define RESP_SIGMA_THRESHOLD 0.13
#define RESP_PEAK_THRESHOLD 0.1

/**
 * Main routine for peak detection on a dataseries
 *
 * @param data      Input dataseries, previosly loaded from file
 * @param n         Length of dataseries
 *
 * The following are optional inputs. These can be given so they can be saved
 * to file once ampd is done. If nullpointer is given, ampd malloc memory for
 * subsequent calculations and frees it up once completed. The input data should
 * already be filtered and smoothed if necessary.
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

int ampdcpu(float *data, int n, struct ampd_param *param,
            struct fmtx *lms, float *gam,float *sig, int *pks){

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
    double dampling_rate = param->sampling_rate;
    double a = param->a;
    double rnd_factor = param->rnd_factor;
    check_ampd_param(param); // check if ampd param struct is filled properly
    /*
     * linear fitting subroutine
     * y = a * x + b, where x is time in seconds
     */
    double x, sumx, sumx2, sumxy, sumy, sumy2, denom;
    x = 0.0; sumx = 0.0; sumx2 = 0.0; sumxy = 0.0; sumy = 0.0; sumy2 = 0.0;
    for(i=0;i<n;i++){
        x = (double) i * 1 / sampling_rate;
        sumx += x;  sumx2 += x*x; sumxy += x*data[i];
        sumy += data[i]; sumy2 += data[i]*data[i];     
    }
    denom = (n * sumx2 - sumx * sumx);
    if(denom == 0){
        // singular matrix
        param->fit_a = 0; param->fit_b = 0; param->fit_r = 0; 
        fprintf(stderr, "Cannot fit line to data, qitting.\n");
        return -1;
    }
    param->fit_a = (n * sumxy - sumx * sumy) / denom;
    param->fit_b = (sumy * sumx2 - sumx * sumxy) / denom;
    param->fit_r = (sumxy - sumx * sumy / n) / 
                    sqrt((sumx2 - sumx * sumx / n ) * sumx2 - sumy * sumy/n);
    /*
     * detrending
     *
     */
    for(i=0; i<n; i++)
        data[i] -= (float)(a * (double)i / sampling_rate + b);

    /*
     * calculating LMS
     *
     */
    l = (int)ceil(n/2)-1; // rows of LMS
    if(null_inputs[0] == 1){ 
        // setup lms struct if nullpoier was given as input
        lms = malloc_fmtx(l, n);
    }
    for(i=0; i<n; i++){
        
    }


    return n_pks;


}
/**
 * Check if param struct is filled properly, exit on error.
 */
void check_ampd_param(struct ampd_param *param){

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

