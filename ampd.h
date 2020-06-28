/*
 * ampd.h
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
#include <string.h>
#include <stdbool.h>
#include <unistd.h>

#define MAX_DATA_LEN 100000
#define MAX_PATH_LEN 256
#define VERBOSE_DEFAULT 0
#define OUTPUT_ALL_DEFAULT 0

#define ALPHA 1 // constant factor

/* Util functions */
/*================*/

void printf_help();
int ceiling(double z);

/* Core functions */
/*================*/

/* linear fit to data */
int linreg(float *data, int n, double *a, double *b);
/* subtract least squares fit*/
int linear_detrend(float *data, int n, double a, double b);
/* malloc for local maxima scalogram matrix*/
int malloc_lms(struct Mtx *mtx, int rows, int cols);
/* row summation of local maxima scalogram*/
int row_sum_LMS(struct Mtx *mtx, float *gamma);

/* general matrix */
struct Mtx {

    int rows;
    int cols;
    float **data;
}
