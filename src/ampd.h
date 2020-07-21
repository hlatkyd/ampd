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
 *
 * main AMPD routine is found in ampdr.c
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
 * Default AMPD parameters. change these if needed
 *
 * sampling rate: data sampling rate in Hz
 * sigma_threshold: sigma is counted as zero below this value meaning
 *                  it's a peak at that index
 * peak_threshold:  peaks should not be closer than this, in seconds
 * overlap:         long data (over minutes, or 10k samples) is processed
 *                  in batches, and batches can overlap to help detect
 *                  peaks by going over them multiple times.
 *                  overlap 0 means no overlap, 0.5 means data is processed
 *                  twice, and so on..
 ********************************************************************
 *
 */
//Default AMPD parameters for data agnostic usage
#define DEF_SAMPLING_RATE 100
#define DEF_SIGMA_THRESHOLD 0.1
#define DEF_PEAK_THRESHOLD 0.05
#define DEF_OVERLAP 0

#define DEF_A 1
#define DEF_RND_FACTOR 1

// Default AMPD parameters for respiration
#define RESP_SAMPLING_RATE 100
#define RESP_SIGMA_THRESHOLD 0.1
#define RESP_PEAK_THRESHOLD 0.1

// Default AMPD parameters for pulsoxymetry
#define PULS_SAMPLING_RATE 100
#define PULS_SIGMA_THRESHOLD 0.50
#define PULS_PEAK_THRESHOLD 0.05
/*
 ********************************************************************
 */

#define MAX_PATH_LEN 1024

// defaults
#define VERBOSE 0
#define SMOOTH_DATA 0
#define OUTPUT_ALL 0
#define OUTPUT_LMS 0     // full and reduced local maxima scalogram
#define OUTPUT_RATE 1   // output peaks per min to file

struct ampd_config{

    /* general io*/

    int smooth_data;
    int output_all;
    int output_lms;
    int output_rate;
    double hpfilt;
    double lpfilt;
    /* data specific*/
    double def_sampling_rate;
    double def_sigma_threshold;
    double def_peak_threshold;
    double def_overlap;

    double resp_sampling_rate;
    double resp_sigma_threshold;
    double resp_peak_threshold;
    double resp_overlap;

    double puls_sampling_rate;
    double puls_sigma_threshold;
    double puls_peak_threshold;
    double puls_overlap;

};


void printf_help();
void printf_data(float *data, int n);

/* set ampd_params from defaults*/
//TODO do this from conf file
void set_ampd_param(struct ampd_param *p, char *type);
void set_ampd_param_cfg(struct ampd_param *p, struct ampd_config *cfg);
/* parse config file*/
void preload_config(char *path, struct ampd_config *conf);
void load_config(char *path, struct ampd_config *conf, char *datatype);
/* saving and loading data data*/
int fetch_data(char *path, float *data, int n, int ind);
int mkpath(char *file_path, mode_t mode);
void save_fmtx(struct fmtx *mtx, char *path);
void save_data(void *data, int n, char *path, char *type);
void save_ampd_param(struct ampd_param *param, char *path);

int count_char(char *path, char cc);
/* extract filename from full path and omitting file extension*/
void extract_raw_filename(char *path, char *filename, int bufsize);

//TODO
/* merge peak indices from subsequent batches*/
int merge_peaks(int *sum_peaks, int sum_n, int *peaks, int n, int ind);
