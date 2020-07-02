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
#include <libgen.h>
#include <stdbool.h>
#include <unistd.h>
#include <math.h>

#define MAX_DATA_LEN 1000000    // testing, load this many points only
#define MAX_PATH_LEN 256
#define DATA_BUF 5000  // number of data points to work on at a time
// defaults
#define VERBOSE 0
#define OUTPUT_ALL 0
#define OUTPUT_LMS 0     // full and reduced local maxima scalogram
#define OUTPUT_VECTORS 1 // sigma, gamma, peaks

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
/* extract filename from full path and omitting file extension*/
void extract_raw_filename(char *path, char *filename, int bufsize);

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
void calc_lms(struct Mtx *lms, float *data);
/* row summation of local maxima scalogram*/
void row_sum_lms(struct Mtx *lms, double *gamma);
/* make rescaled LMS */
int rescale_lms(struct Mtx *lms, struct Mtx *rlms, double *gamma, double *gamma_min);
/* calculate column-wise standard dev of rescaled LMS*/
int col_stddev_lms(struct Mtx *rlms, double *sigma, int gamma_min);
/* find peaks: where sigma = 0 */
int find_peaks(double *sigma, int n, int *peaks, int *n_peaks);
/* main ampd routine*/
int ampd(float *data, int n, struct Mtx *lms, struct Mtx *rlms,
         double *gamma, double *sigma, int *peaks);
/* concatenate peak indices from subsequent batches*/
int concat_peaks(int *sum_peaks, int *peaks, int ind);
