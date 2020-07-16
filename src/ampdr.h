/*
 * ampdr.c
 *
 * Implementation of AMPD for CPU with some optimization .
 *
 */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

/*
 * Default AMPD parameters for data agnostic usage
 */
#define DEF_SAMPLING_RATE 200
#define DEF_SIGMA_THRESHOLD 0.1
#define DEF_PEAK_THRESHOLD 0.1
#define DEF_OVERLAP 0
#define DEF_
/*
 * Default AMPD parameters for respiration
 */
#define RESP_SAMPLING_RATE 200
#define RESP_SIGMA_THRESHOLD 0.1
#define RESP_PEAK_THRESHOLD 0.1
/*
 * Default AMPD parameters for pulsoxymetry
 */
#define PULSOX_SAMPLING_RATE 200
#define PULSOX_SIGMA_THRESHOLD 0.1
#define PULSOX_PEAK_THRESHOLD 0.1



// for long only inputs
#define OUTPUT_ALL 10
#define OUTPUT_LMS 20
#define OUTPUT_RATE 30
#define OVERLAP 40

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
void set_ampd_param(struct ampd_param *p, char *type);

