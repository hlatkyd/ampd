/*
 * ampdr.c
 *
 * Implementation of AMPD for CPU with some optimization .
 *
 */

#include "ampdr.h"

/**
 * Main routine for peak detection on a dataseries
 *
 * @param data      Input dataseries, previosly loaded from file
 * @param n         Length of dataseries
 * @param param     pointer to struct of AMPD parameters, such as sampling rate,
 *                  thresholds, etc
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
            struct fmtx *lms, double *gamma,double *sigma, int *pks){

    /* To keep track of nullpointer inputs. 1 means nullponter input. in order:
     * lms, gam, sig, pks
     */
    int null_inputs[4] = {0,0,0,0}; 
    if(lms == NULL)
        null_inputs[0] = 1;
    if(gamma == NULL)
        null_inputs[1] = 1;
    if(sigma == NULL)
        null_inputs[2] = 1;
    if(pks == NULL)
        null_inputs[3] = 1;
    int i, k;
    double sampling_rate = param->sampling_rate;
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
    param->fit_r = (sumxy - sumx * sumy / n) / \
                    sqrt((sumx2 - sumx * sumx / n ) * sumx2 - sumy * sumy/n);
    /*
     * detrending
     *
     */
    for(i=0; i<n; i++)
        data[i]-=(float)(param->fit_a*(double)i / sampling_rate+param->fit_b);

    /*
     * calculating LMS
     *
     */
    int l = (int)ceil(n/2)-1; // rows of LMS
    float rnd;
    double rnd_factor = param->rnd_factor;
    double a = param->a;
    if(null_inputs[0] == 1){ 
        // setup lms struct if nullpointer was given as input
        lms = malloc_fmtx(l, n);
    }
    for(i=0; i<n; i++){
        for(k=0; k<l; k++){
            rnd = (float) rand() / (float)RAND_MAX * (float)rnd_factor;
            if(i<k || i>n-k+1){
                lms->data[k][i] = rnd + a;
                continue;
            }
            if(data[i-1] > data[i-k-1] && data[i-1] > data[i+k-1])
                lms->data[k][i] = 0.0;
            else
                lms->data[k][i] = rnd + a;
        }
    }
    /*
     * calculating gamma, and find its minimum, lambda
     */
    double min = 0.0;
    int lambda = 0;
    if(null_inputs[1] == 1){
        gamma = malloc(sizeof(double) * l);
    }
    for(k=0; k<l; k++){
        if(k==1)
            min = gamma[0]; // set first
        gamma[k] = 0.0;
        for(i=0; i<n; i++){
            gamma[k] += lms->data[k][i];
        }    
        if(gamma[k] < min){
            min = gamma[k];
            lambda = k;
        }
    }
    param->lambda = lambda;
    /*
     * calculating sigma and find the peaks
     */
    double sum_m_i;
    double sum_outer;
    double sigma_thresh = param->sigma_thresh;
    int ind_thresh = param->peak_thresh / param->sampling_rate;
    int n_pks = 0; int j=0;
    if(null_inputs[2] == 1)
        sigma = malloc(sizeof(double) * n);
    if(null_inputs[3] == 1)
        pks = malloc(sizeof(int)*n);

    for(i=0; i<n; i++){
        sigma[i] = 0.0;
        sum_m_i = 0.0;
        for(k=0; k<lambda; k++)
            sum_m_i += lms->data[k][i] / (double) lambda;
        for(k=0; k<lambda; k++)
            sigma[i] += fabs(lms->data[k][i]-sum_m_i) / (double)(lambda-1);
        // check for peak
        if(sigma[i] < sigma_thresh){
            if(i - pks[j-1] > ind_thresh){
                pks[j] = i;
                j++;
                n_pks++;
            }
            else
                continue;
        }
    }
    // free memory if aux output is not needed
    if(null_inputs[0] == 1){
        for(i=0; i<l; i++)
            free(lms->data[j]);
        free(lms);
    }
    if(null_inputs[1] == 1)
        free(gamma);
    if(null_inputs[2] == 1)
        free(sigma);
    if(null_inputs[3] == 1)
        free(pks);
    return n_pks;
}

struct fmtx *malloc_fmtx(int rows, int cols){

    struct fmtx *mtx = malloc(sizeof(struct fmtx));
    mtx->rows = rows;
    mtx->cols = cols;
    mtx->data = malloc((mtx->rows*sizeof(float *)));
    for(int i=0; i<mtx->rows;i++){
        mtx->data[i] = malloc(mtx->cols * sizeof(float));
    }
    return mtx;

}

