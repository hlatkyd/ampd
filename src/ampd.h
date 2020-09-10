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
#define VERSION_MAJ 0
#define VERSION_MIN 5

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
 * Default AMPD parameters.
 * Fallback defaults, in case these are not specified at input and
 * no config file is used.
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
#define DEF_BATCH_LENGTH 60

#define DEF_SIGMA_THRESHOLD 0.01
#define DEF_PEAK_THRESHOLD 0.05

#define DEF_A 1             // works ok, do not change
#define DEF_RND_FACTOR 1    // works ok, do not change
#define DEF_N_ZPAD 50       // TODO, padding,currently unused
#define DEF_N_BINS 50       // histogram bins for data flipping

// Default AMPD parameters for respiration
#define RESP_SAMPLING_RATE 100
#define RESP_SIGMA_THRESHOLD 0.01
#define RESP_PEAK_THRESHOLD 0.1

// Default AMPD parameters for pulsoxymetry
#define PULS_SAMPLING_RATE 100
#define PULS_SIGMA_THRESHOLD 0.03
#define PULS_PEAK_THRESHOLD 0.05
/***************************************************************************/
/* Default preprocess parameters.
 * Preprocess does  the same on the batches as the utility program ampdpreproc
 * does on the entire data.
 *
 */
#define RESP_PREPROC 1
#define RESP_HPFILT 0.5
#define RESP_LPFILT 3

#define PULS_PREPROC 1
#define PULS_HPFILT 2
#define PULS_LPFILT -1

#define DEF_PREPROC 1
#define DEF_HPFILT 0.1
#define DEF_LPFILT -1

/*************************************************************************/

// These are booleans: set either 0 or 1
// output peaks per minutes to file for each window
#define DEF_OUTPUT_RATE 1
// output peak indices to file
#define DEF_OUTPUT_PEAKS 1
// output metadata to file
#define DEF_OUTPUT_META 1
// output aux files, useful for troubleshooting
#define DEF_OUTPUT_ALL 0
// output local maxima scalogram
#define DEF_OUTPUT_LMS 0
// TODO output plot images from aux data
#define DEF_OUTPUT_IMG 0

#define MAX_PATH_LEN 1024


// only used for saving to meta file
struct meta_param{

    char infile[MAX_PATH_LEN];
    char basename[MAX_PATH_LEN];
    char datatype[32];
    double sampling_rate;
    int preproc;
    double hpfilt;
    double lpfilt;
    double batch_length;
    int total_batches;
    int total_peaks;
};

// settings for preprocessing: smooothing and filtering
struct preproc_param{

    int preproc;
    double lpfilt;
    double hpfilt;

};

struct batch_param{

    int ind;    // index of current batch
    int cycles;// total number of batches
    int n; // same as data array length
    int l;
    double batch_length; // data length in seconds
    double sampling_rate;   // sampling rate in Hz
    int n_peaks;    // peak count in batch
    double peaks_per_min;
    double mean_pk_dist;
    double stdev_pk_dist;


};
//TODO
// this is pretty much unused, cleaup or finish needed
struct ampd_config{

    /* general io*/
    int smooth_data;
    int output_all;
    int output_lms;
    int output_rate;
    int output_peaks;

    /* data specific*/
    int preprocess;
    char datatype[32];
    double hpfilt;
    double lpfilt;

    double sampling_rate;
    double batch_length;
    double sigma_threshold;
    double peak_threshold;
    double overlap;

};


void printf_help();
void printf_data(float *data, int n);

/* initialize parameter structs*/
void init_preproc_param(struct preproc_param *p);
void init_ampd_config(struct ampd_config *cfg);
/* set ampd_params from defaults*/

void set_preproc_param(struct preproc_param *p, char *type);
int fork_ampdpreproc(char *path, double lpfilt, double hpfilt);
void set_ampd_param(struct ampd_param *p, char *type);
/* parse config file*/
void preload_config(char *path, struct ampd_config *conf);
void load_config(char *path, struct ampd_config *conf, char *datatype);
/* saving and loading data data*/
int fetch_data(char *path, float *data, int n, int ind, int n_zpad);
/* get data batches from full data array*/
int fetch_data_buff(float *full_data, int len,float *data,int n,int ind,int n_zpad);
/* load data from file to memory*/
void load_from_file(char *path, float *full_data, int n);

/* make histogram to flip the data in case inhales are minima*/
void histogram(float *data, int n, int *bins, int n_bins);
double centre_of_mass(int *bins, int n_bins);
void flip_data(float *data, int n);

int mkpath(char *file_path, mode_t mode);
void save_fmtx(struct fmtx *mtx, char *path);
void save_data(void *data, int n, char *path, char *type);
void save_rate(double rate, char *path);
void save_ampd_param(struct ampd_param *param, char *path);
void save_batch_param(struct batch_param *p, char *path);
void save_meta(struct meta_param *p, struct preproc_param *pp, char *path);

int count_char(char *path, char cc);
/* extract filename from full path and omitting file extension*/
void extract_raw_filename(char *path, char *filename, int bufsize);

// UNUSED
/*-----------------------------------------------------------------*/
//TODO
/* merge peak indices from subsequent batches*/
int merge_peaks(int *sum_peaks, int sum_n, int *peaks, int n, int ind);

//TODO do this from conf file
void set_ampd_param_cfg(struct ampd_param *p, struct ampd_config *cfg);
