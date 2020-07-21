/*
 * ampdr.c
 *
 * Implementation of AMPD for CPU with some optimization .
 *
 * AMPD
 * ====
 * Peak detection algorithm for quasiperiodic data. Main usage of this
 * implementation is detection of peaks in rat physiological data:
 * respiration and pulsoxymmetry waveforms.
 *
 * Reference paper:
 * An Efficient Algorithm for Automatic Peak Detection in Noisy
 * Periodic and Quasi-Periodic Signals
 * DOI:10.3390/a5040588

 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

/* generic matrix of float */
struct fmtx {

    int rows;
    int cols;
    float **data;

};

struct ampd_param {

    double sampling_rate;
    char datatype[128];
    /* ampd constant factors for LMS calculation*/
    double a;               // alpha as in reference paper
    double rnd_factor;      // multiplier of rand[0,1]
    /* linear fitting result params y = a*x + b; r is residual*/
    double fit_a;
    double fit_b;
    double fit_r;
    int lambda;          // reduced LMS lambda
    double sigma_thresh;    // sigma threshold above 0
    double peak_thresh;     // peak minimum distance in seconds

};
/* main routine */
int ampdcpu(float *data,int n, struct ampd_param *param, 
            struct fmtx *lms,double *gam, double *sig, int *pks);

/* helper routines */
int linregu(float *y, int n, double rate, double *a, double *b, double *r);
/* find lambda*/
int more_sophisticated_way_to_lambda(double *gamma, int l);

/* util */
struct fmtx *malloc_fmtx(int rows, int cols);

