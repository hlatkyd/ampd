/*
 * ampdr.c
 *
 * Implementation of AMPD for CPU with some optimization .
 *
 */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define DATATYPE RESP

/* Respiration peak counting constants*/
#define RESP 0
#define RESP_TIMERES 0.0002
/* Pulsoxy peak counting constants*/
#define PULSOX_TIMERES 0.0002
#define PULSOX 1

#define ECG 2

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


/* util */
struct fmtx *malloc_fmtx(int rows, int cols);
int set_param(struct ampd_param *param);

