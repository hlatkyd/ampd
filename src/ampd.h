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
#include <getopt.h>
#include <errno.h>
#include <assert.h>
#include <unistd.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <math.h>
#include <time.h>

#include "ampdr.h"

/*
 * Important parameters, change these to fine-tune for an application
 *
 * TOLERANCE:
 *
 * A tolerance value for peak detection, as sigma is not
 * exacly zero. Some usual values:
 * Rat respiration: 0.10
 *
 * PEAK_MIN_DIST:
 *
 * Hard threshold to ignore assumed peaks which are too close to each other.
 * It is in units of seconds. Usual values:
 * Rat respiration: 0.1
 *
 * SMOOTH_TIMECONST
 *
 * Data is smoothed with moving average. This is the width of the window in 
 * seconds.
 *
 */
#define TOLERANCE 0.13
#define PEAK_MIN_DIST 0.1       // not used, threshold by time distance
#define IND_THRESH 5            // hard threshold by indice distance
#define SMOOTH_TIMEWINDOW 0.005
#define TIMESTEP_DEFAULT 0.01 // sampling time in sec

#define MAX_PATH_LEN 1024

// defaults
#define VERBOSE 0
#define SMOOTH_DATA 1
#define OVERLAP_DEF 0.0 // overlapping batches in time domain
#define DATA_BUF_DEF 5000  // number of data points to work on at a time
#define OUTPUT_ALL 0
#define OUTPUT_LMS 0     // full and reduced local maxima scalogram
#define OUTPUT_RATE 1   // output peaks per min to file
#define OUTPUT_VECTORS 1 // sigma, gamma, peaks

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

int mkpath(char *file_path, mode_t mode);
void save_mtx(struct Mtx *mtx, char *path);
void save_data(float *data, int n, char *path);
void save_ddata(double *data, int n, char *path);
void save_idata(int *data, int n, char *path);
void save_fitdata(double a, double b, double n, char *path);

int count_char(char *path, char cc);
/* extract filename from full path and omitting file extension*/
void extract_raw_filename(char *path, char *filename, int bufsize);
/* return the index of the global minumum of a vector*/
int argmin_minind(double *data, int n);
int calc_halfwindow(double timestep, double timewindow);

/* Core functions */
/*================*/

/* load a paort of data from a file*/
int fetch_data(char *path, float *data, int n, int ind );
/* smooth data */
void smooth_data(float *data, int n, int wh, float *newdata, int new_n);
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
int rescale_lms(struct Mtx *lms,struct Mtx *rlms,double *gamma,double *gamma_min);
/* calculate column-wise standard dev of rescaled LMS*/
void col_stddev_lms(struct Mtx *rlms, double *sigma, int gamma_min);
/* find peaks: where sigma = 0 */
void find_peaks(double *sigma, int n, int *peaks, int *n_peaks);

void catch_false_peaks(int *peaks, int *n_peaks, double timestep, double thresh);
/* main ampd routine*/
int ampd(float *data,int n,struct Mtx *lms,double *gamma,double *sigma,int *peaks,double *a, double *b);

/* concatenate peak indices from subsequent batches*/
int concat_peaks(int *sum_peaks, int sum_n, int *peaks, int n, int ind);

/* TODO
 * An more integrated, optimized version of ampd.
 * Returns number of peaks.
 */
int ampd2(float *data,int n,struct Mtx *lms,double *gam,double *sig,int *pks);
