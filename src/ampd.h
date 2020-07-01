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
#include <stdlib.h>
#include <string.h>
#include <stdbool.h>
#include <unistd.h>
#include <math.h>

#define MAX_DATA_LEN 1000000    // testing, load this many points only
#define MAX_PATH_LEN 256
#define DATA_BUF 50  // number of data points to work on at a time
#define VERBOSE_DEFAULT 0
#define OUTPUT_ALL_DEFAULT 0

// time between samples in seconds
// for wighing, to avoid overflow effects
#define TIMESTEP_DEFAULT 0.0001 

#define ALPHA 1 // constant factor

/*****************
 * general matrix
 *
 * indexing: Mtx[row][col]
 */

struct Mtx {

    int rows;
    int cols;
    float **data;
};

/* Util functions */
/*================*/

void printf_help();
void printf_data(float *data, int n);
void save_mtx(struct Mtx *mtx, char *path);
void save_data(float *data, int n, char *path);
int count_char(char *path, char cc);

/* Core functions */
/*================*/

int fetch_data(char *path, float *data, int n, int ind);
/* linear fit to data */
int linear_fit(float *data, int n, double ts, double *a, double *b, double *r);
/* subtract least squares fit*/
int linear_detrend(float *data, int n, double ts, double a, double b);
/* malloc for a general matrix*/
struct Mtx* malloc_mtx(int rows, int cols);
/* calculate local maxima scalogram */
int calc_lms(struct Mtx *lms, float *data);
/* row summation of local maxima scalogram*/
int row_sum_LMS(struct Mtx *mtx, float *gamma);

