/*
 * ampdr.c
 *
 * Implementation of AMPD for CPU with some optimization .
 *
 */


#define DATATYPE RESP

/* Respiration peak counting constants*/
#define RESP 0
#define RESP_TIMERES 0.0002
/* Pulsoxy peak counting constants*/
#define PULSOX_TIMERES 0.0002
#define PULSOX 1

#define ECG 2
/* matrix of float */
struct fmtx {
    int rows;
    int cols;
    float **data;
};
/* main routine */
int ampdcpu(float *data,int n,struct fmtx *lms,float *gam,float *sig, int *pks);

void set_constants(char *type, double *tol, int *min_dst);
