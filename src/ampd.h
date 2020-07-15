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
#include "filters.h"

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
#define PEAK_MIN_DIST 0.05       // not used, threshold by time distance
#define IND_THRESH 10            // hard threshold by indice distance
#define SMOOTH_TIMEWINDOW 0.005
#define TIMESTEP_DEFAULT 0.01 // sampling time in sec

#define MAX_PATH_LEN 1024

// defaults
#define VERBOSE 0
#define SMOOTH_DATA 0
#define OVERLAP_DEF 0.0 // overlapping batches in time domain
#define DATA_BUF_DEF 5000  // number of data points to work on at a time
#define OUTPUT_ALL 0
#define OUTPUT_LMS 0     // full and reduced local maxima scalogram
#define OUTPUT_RATE 1   // output peaks per min to file
#define OUTPUT_VECTORS 1 // sigma, gamma, peaks

#define ALPHA 1 // constant factor
#define RAND_FACTOR 1


/* Util functions */
/*================*/

void printf_help();
void printf_data(float *data, int n);

void set_ampd_param(struct ampd_param *p);
/* saving and loading data data*/
int fetch_data(char *path, float *data, int n, int ind);
int mkpath(char *file_path, mode_t mode);
void save_fmtx(struct fmtx *mtx, char *path);
void save_data(void *data, int n, char *path, char *type);

int count_char(char *path, char cc);
/* extract filename from full path and omitting file extension*/
void extract_raw_filename(char *path, char *filename, int bufsize);

/* merge peak indices from subsequent batches*/
int merge_peaks(int *sum_peaks, int sum_n, int *peaks, int n, int ind);
