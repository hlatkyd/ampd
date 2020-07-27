/*
 * ampd.c
 *
 */

#include "ampd.h"

#define TESTING 0
#define CONF_PATH "ampd.conf"

/* only for getopt */
#define ARG_OUTPUT_ALL 1
#define ARG_OUTPUT_LMS 2
#define ARG_OUTPUT_RATE 3
#define ARG_OUTPUT_PEAKS 8
#define ARG_OUTPUT_IMG 9
#define ARG_OVERLAP 4
#define ARG_PREPROC 5
#define ARG_HPFILT 6
#define ARG_LPFILT 7

int verbose;
int output_all = DEF_OUTPUT_ALL;   // output all intermediary data except LMS
int output_lms = DEF_OUTPUT_LMS;   // output local maxima scalogram, HUGE!
int output_rate = DEF_OUTPUT_RATE;  // output peaks per min
int output_peaks = DEF_OUTPUT_PEAKS; // output peak indices
int output_img = DEF_OUTPUT_IMG; // save plot image in all batches for inspection
int preproc = DEF_PREPROC; 
static struct option long_options[] = 
{
    {"infile",required_argument, NULL, 'f'},
    {"outdir",required_argument, NULL, 'o'},
    {"auxdir",required_argument, NULL, 'a'},
    {"verbose", optional_argument, NULL, 'v'},
    {"help",no_argument, NULL, 'h'},
    {"datatype",required_argument,NULL, 't'},
    {"sampling-rate",required_argument, NULL, 'r'},
    {"batch-length", required_argument, NULL, 'l'},
    {"overlap", required_argument, NULL, ARG_OVERLAP},
    /* preprocess with ampdpreproc*/
    {"preproc", no_argument, NULL, ARG_PREPROC},
    {"hpfilt", required_argument, NULL, ARG_HPFILT},
    {"lpfilt", required_argument, NULL, ARG_LPFILT},
    {"output-all", no_argument, NULL, ARG_OUTPUT_ALL},
    {"output-lms", no_argument, NULL, ARG_OUTPUT_LMS},
    {"output-rate", no_argument, NULL, ARG_OUTPUT_RATE}, // unused
    {"output-peaks", no_argument, NULL, ARG_OUTPUT_PEAKS},
    {"output-img", no_argument, NULL, ARG_OUTPUT_IMG}, //TODO
    {NULL, 0, NULL, 0}
};
static char optstring[] = "hvf:o:a:l:r:t:";

/**
 * Print description and general usage
 */
void printf_help(){

    printf(
    "ampd - Automatic Multi-scale Peak detection\tv%d.%d\n"
    "=========================================================================\n"
    ,VERSION_MAJ,VERSION_MIN);
    /*
    "Peak detection algorithm for quasiperiodic data. Main usage of this "
    "implementation is detection of peaks in rat physiological data: "
    "respiration and pulsoxymmetry waveforms.\n\n"
    */
    /*
    "Reference paper:\n"
    "An Efficient Algorithm for Automatic Peak Detection in Noisy "
    "Periodic and Quasi-Periodic Signals\n"
    "DOI:10.3390/a5040588\n\n"
    */
    printf(
    "This program takes a file input which contains a single timeseries "
    "of quasi-periodic data. The output is a file containing the indices "
    "of the peaks as calculated.\n"

    "Input file should only contain a float value in each line.\n"
    "Main output files contains the indices of peaks, while aux output "
    "directory contains various intermediate data for error checking. "
    "The final peak count is sent to stdout."
    "\n"

    "\n"
    // TODO cleanup
    "Usage from command line:\n\n"
    "ampd -f [infile] -o [out_dir] -r [sampling_rate] -l [batch_length]\n\n"
    "Optional arguments:\n"
    "-o --outdir:\t\toutput dir\n"
    "-v --verbose:\t\tverbose\n"
    "-h --help:\t\tprint help\n"
    "-r --samplig-rate:\tsampling rate input data in Hz\n"
    "-l --batch-length:\tdata window length in seconds\n"
    "--preproc:\t\tcall ampdpreproc on data first\n"
    "--lpfilt:\t\tapply lowpass filter in Hz\n"
    "--hpfilt:\t\tapply highpass filter in Hz\n"
    "-a --auxdir:\t\taux data root dir, default is cwd/ampd.aux\n"
    //"--overlap:\tmake batches overlapping in time domain\n"
    "--output-all:\t\toutput aux data, local maxima scalogram not included\n"
    "--output-lms:\t\toutput local maxima scalogram (high disk space usage)\n"
    "--output-rate:\t\toutput peak-per-min\n"
    "--output-peaks\t\toutput peak indices\n"
    "\n"
            );
}
/**
 * Main handles the command line input parsing and output file management.
 */
int main(int argc, char **argv){

    int opt;
    struct ampd_config *conf;
    char datatype[32] = {0};
    // timing
    clock_t begin, end;
    double time_spent;

    char infile[MAX_PATH_LEN] = {0};
    char infile_basename[MAX_PATH_LEN]; // input, without dir and extension
    char outdir[MAX_PATH_LEN] = {0};
    char outfile_peaks[MAX_PATH_LEN] = {0}; // main output file with indices of peaks
    char outfile_rate[MAX_PATH_LEN] = {0};
    FILE *fp_out;        // main output file containing the peak indices
    FILE *fp_out_rate;
    char cwd[MAX_PATH_LEN]; // current directory
    /*
     * aux output files
     */
    char raw_path[MAX_PATH_LEN];
    char detrend_path[MAX_PATH_LEN];
    char lms_path[MAX_PATH_LEN];
    char rlms_path[MAX_PATH_LEN];
    char gamma_path[MAX_PATH_LEN];
    char sigma_path[MAX_PATH_LEN];
    char peaks_path[MAX_PATH_LEN];
    char preproc_path[MAX_PATH_LEN];
    char param_path[MAX_PATH_LEN];
    char batch_param_path[MAX_PATH_LEN];
    char rate_path[MAX_PATH_LEN];

    /* load data for preproc*/
    float *full_data;
    /*
     * batch processing
     */
    float *data;        //  timeseries
    int n;              // number of elements in timeseries, in a batch, dynamic
    int ind;            // index of next batch
    int n_zpad = DEF_N_ZPAD; // zeropad data to help detect peaks on edges
    int data_buf;
    int datalen;        // full data length
    double batch_length = DEF_BATCH_LENGTH;
    double overlap = 0.0;          //overlap ratio
    int n_overlap;             // overlapping data points
    int cycles;         // number of data batches
    int sum_n_peaks;    // summed peak number from all batches
    int *sum_peaks;      // concatenated vector from peak indices
    double peaks_per_min;
    double data_buf_sec = 0.0;
    struct batch_param *bparam; // only for outputting batch utility parameters

    /* 
     * filtering
     */
    struct preproc_param *pparam;
    /*
     * ampd routine
     */
    struct ampd_param *param;
    double sampling_rate = 0.0;
    int l;
    struct fmtx *lms;
    double *gamma;
    double *sigma;
    int *peaks;
    int n_peaks;        // main output, number of peaks

    // main output file base
    char outdir_def[] = "ampd.out"; //
    // peaks will be save to ampd.out.pks
    // peaks per min to ampd.out.rate
    char outbase[256];

    // aux output paths
    char aux_dir[MAX_PATH_LEN] = {0};
    char batch_dir[MAX_PATH_LEN] = {0};
    char aux_dir_def[] = "ampd.aux"; // full default is cwd plus this
    if(argc == 1){
        printf_help();
        return 0;
    }
    // load config first, inputs can overwrite these
    char conf_path[] = CONF_PATH;
    conf = malloc(sizeof(struct ampd_config));
    param = malloc(sizeof(struct ampd_param));
    bparam = malloc(sizeof(struct batch_param));
    pparam = malloc(sizeof(struct preproc_param)); // filtering params
    memset(bparam, 0, sizeof(bparam));
    memset(conf, 0, sizeof(struct ampd_config));
    init_preproc_param(pparam);

    set_ampd_param(param, "def"); // ampd_param is for the AMPD routine only
    //TODO
    //load_config(conf_path, conf, NULL);

    // parse options
    while((opt = getopt_long(argc,argv,optstring,long_options,NULL)) != -1){
        switch(opt){
            case 'h':
                printf_help();
                return 0;
            case 'v':
                if(optarg != NULL){
                    verbose = atoi(optarg);
                }
                else {
                    verbose = 1;
                }
                break;
            case 't':
                strcpy(datatype,optarg);
                strcpy(conf->datatype, optarg);
                break;
            case 'o':
                strcpy(outdir, optarg);
                break;
            case 'f':
                strcpy(infile, optarg);
                break;
            case 'r':
                sampling_rate = atof(optarg);
                conf->sampling_rate = sampling_rate;
                break;
            case 'a':
                strcpy(aux_dir, optarg);
                break;
            case 'l':
                batch_length = atof(optarg);
                conf->batch_length = batch_length;
                break;
            case ARG_OUTPUT_ALL:
                output_all = 1;
                conf->output_all = 1;
                break;
            case ARG_OUTPUT_LMS:
                output_lms = 1;
                conf->output_lms = 1;
                break;
            case ARG_OUTPUT_PEAKS:
                output_peaks = 1;
                conf->output_peaks = 1;
                break;
            case ARG_OUTPUT_RATE:
                output_rate = 1;
                conf->output_rate = 1;
                break;
            case ARG_OVERLAP:
                overlap = atof(optarg);
                if(overlap > 1){
                    printf("error: overlap higher than 1.0 not permitted\n");
                    exit(EXIT_FAILURE);
                }
                //TODO remove once done
                fprintf(stderr,"Warning: overlap > 0 is not finished yet.\n");
                fprintf(stderr,"Defaulting to  overlap = 0.\n");
                overlap = 0;
                break;
            case ARG_PREPROC:
                // do filtering and smoothing
                preproc = 1;
                pparam->preproc = 1;
                break;
            case ARG_HPFILT:
                pparam->hpfilt = atof(optarg);
                break;
            case ARG_LPFILT:
                pparam->lpfilt = atof(optarg);
                break;
            case ARG_OUTPUT_IMG:
                output_img = 1;
                break;
            case '?':
                break;
            case ':':
                break;

        }
    }
    /* Setting up output paths.
     * Main output path is given as command line argument, if not it is the cwd
     * General format of aux output files:
     * aux_root/data_ID/batch_[i]/out_intermediate_files
     * aux_root is the current working dir by default, but can be set with a
     * command line argument. data_ID contains the name of the file where the
     * input data is coming from.
     */
    begin = clock();
    getcwd(cwd, sizeof(cwd));

    set_preproc_param(pparam, datatype);
    if(sampling_rate != 0)
        param->sampling_rate = sampling_rate;
    extract_raw_filename(infile, infile_basename, sizeof(infile_basename));
    // setting outptu files
    snprintf(outfile_peaks, sizeof(outfile_peaks), "%s/%s/%s.peaks",
                cwd, outdir_def,infile_basename);
    snprintf(outfile_rate, sizeof(outfile_rate), "%s/%s/%s.rate",
                cwd, outdir_def,infile_basename);
    if(strcmp(aux_dir,"")==0){
        snprintf(aux_dir, sizeof(aux_dir), "%s/%s",cwd, aux_dir_def);
    }
    if(data_buf_sec != 0.0){
        printf("here!\n");
        data_buf = (int)(data_buf_sec * sampling_rate);
    }
    // setting available param
    if(batch_length != 0)
        bparam->batch_length =  batch_length;
    /*
     * set available config
     */
    if(strcmp(datatype,"")==0){
        strcpy(datatype,"def");
        strcpy(conf->datatype,datatype);
    }
    // setting remaining variables for processing
    sum_n_peaks = 0;
    datalen = count_char(infile, '\n');
    data_buf = (int)(batch_length * param->sampling_rate);
    n_overlap = (int)((double)data_buf * overlap);
    cycles = (int) (ceil(datalen / (double) (data_buf - n_overlap)));
    // preload data and apply filters
    full_data = malloc(sizeof(float) * datalen);
    load_from_file(infile, full_data, datalen);
    if(preproc == 1){
        if(pparam->hpfilt > 0){
            tdhpfilt(full_data, datalen, sampling_rate, pparam->hpfilt);
        }
        if(pparam->lpfilt > 0){
            tdlpfilt(full_data, datalen, sampling_rate, pparam->lpfilt);
        }
    }
    // make path
    // opening main output files
    if(output_peaks == 1){
        mkpath(outfile_peaks, 0777);
        fp_out = fopen(outfile_peaks, "w");
        if(fp_out == NULL){
            fprintf(stderr, "cannot open file for writing %s\n",outfile_peaks);
            exit(EXIT_FAILURE);
        }
    }
    if(output_rate == 1){
        mkpath(outfile_rate, 0777);
        fp_out_rate = fopen(outfile_rate, "w");
        if(fp_out_rate == NULL){
            fprintf(stderr,"cannot open file %s\n",outfile_rate);
            exit(EXIT_FAILURE);
        }

    }
    if(verbose > 0){
        printf("ampd input\n-----------------\n");
        printf("verbose: %d\n",verbose);
        printf("infile: %s\n", infile);
        printf("outfile: %s\n", outfile_peaks);
        printf("datatype: %s\n",datatype);
        printf("lowpassfilt=%lf\n",pparam->lpfilt);
        printf("highpassfilt=%lf\n",pparam->hpfilt);
        printf("aux_dir: %s\n", aux_dir);
        printf("sampling_rate: %.5lf\n",param->sampling_rate);
        printf("batch_length: %lf\n",batch_length);
        printf("overlap: %lf, n_overlap: %d\n",overlap, n_overlap);
        printf("datalen: %d\n", datalen);
        printf("data_buf: %d\n",data_buf);
        printf("cycles: %d\n", cycles);
        printf("output-all: %d\n",output_all);
        printf("output-lms: %d\n",output_lms);
        printf("output-rate: %d\n",output_rate);

    }
    /* fill batch param */
    bparam->cycles = cycles;
    bparam->batch_length = batch_length;
    bparam->sampling_rate = sampling_rate;
    /*
     * Processing
     *
     */
    for(int i=0; i<cycles; i++){
        if(TESTING == 1){
            if(i>0)
                break;
        }

        // settings aux output paths
        snprintf(batch_dir, sizeof(batch_dir),"%s/batch_%d",aux_dir,i);
        snprintf(detrend_path,sizeof(detrend_path),"%s/detrend.dat",batch_dir);
        snprintf(raw_path,sizeof(raw_path),"%s/raw.dat",batch_dir);
        snprintf(lms_path,sizeof(lms_path),"%s/lms.dat",batch_dir);
        snprintf(rlms_path,sizeof(lms_path),"%s/rlms.dat",batch_dir);
        snprintf(gamma_path, sizeof(gamma_path),"%s/gamma.dat",batch_dir);
        snprintf(sigma_path, sizeof(sigma_path),"%s/sigma.dat",batch_dir);
        snprintf(peaks_path, sizeof(peaks_path),"%s/peaks.dat",batch_dir);
        snprintf(preproc_path, sizeof(preproc_path),"%s/smoothed.dat",batch_dir);
        snprintf(param_path, sizeof(param_path),"%s/param.txt",batch_dir);
        snprintf(batch_param_path, sizeof(batch_param_path),"%s/bparam.txt",batch_dir);

        ind = i * (int)(data_buf - n_overlap);
        if(i == cycles-1){
            break; //TODO del test
            //n_init = datalen - i * (int)data_buf;
            n = datalen - i * (int)data_buf;
        }
        else {
            //n_init = (int)data_buf;
            n = (int)data_buf;
        }
        // init matrix and other array pointers
        l = (int)ceil(n/2)-1;
        lms = malloc_fmtx(l, n);
        gamma = malloc(sizeof(double)*l);
        sigma = malloc(sizeof(double)*n);
        peaks = malloc(sizeof(int)*n);

        // fill some more batch param
        bparam->ind = i;
        bparam->n = n;
        bparam->l = l;
        // getting data
        data = malloc(sizeof(float)*n);
        if(verbose > 1){
            printf("\nfetchig data:\n");
            printf("infile=%s\n",infile);
            printf("data=%p\n",data);
            printf("n=%d, ind=%d\n",n,ind);
        }
        //TODO clean this up
        //fetch_data(infile, data, n, ind, n_zpad);
        fetch_data_buff(full_data, data, n, ind, n_zpad);
        if(output_all == 1)
            save_data(data, n, raw_path,"float"); // save raw data

        /* Filter data depending on datatype. Pulsoxy signals benefit from high
         * pass filtering near the respiration frequency.
         *
         */
        /*
        if(pparam->preproc == 1){
            if(pparam->lpfilt > 0.0){
                tdlpfilt(data, n, sampling_rate, pparam->lpfilt);
            }
            if(pparam->hpfilt > 0.0){
                tdhpfilt(data, n, sampling_rate, pparam->hpfilt);
            }
        }
        */
        if(output_all == 1)
            save_data(data, n, preproc_path,"float"); // save smoothed data

        // main ampd routine, contains malloc
        n_peaks = ampdcpu(data, n, param, lms, gamma, sigma, peaks);

        //catch_false_pks(peaks, &n_peaks, ts, thresh);
        sum_n_peaks += (int)(n_peaks * (1.0-overlap));
        if(verbose > 0){
            printf("batch=%d/%d, n_peaks=%d, sum_n_peaks=%d\n",
                    i,cycles, n_peaks,sum_n_peaks);
        }
        // calc peak rate
        bparam->n_peaks = n_peaks;
        peaks_per_min = (double)n_peaks / batch_length * 60.0;
        bparam->peaks_per_min = peaks_per_min; 

        //concat_peaks(sum_peaks, peaks, ind);
        //TODO
        // save aux
        //

        if(output_rate == 1){
            fprintf(fp_out_rate,"%d\n",(int)peaks_per_min);
        }
        if(output_peaks == 1){
            //save_peaks(peaks_path);
            for(int j=0;j<n_peaks;j++){
                fprintf(fp_out,"%d\n",peaks[j]+ind);
            }
        }
        if(output_all == 1){
            //save_data(data, n, raw_path); // save raw data
            save_data(data, n, detrend_path,"float"); // save detrended data
            save_data(sigma, n, sigma_path, "double");
            save_data(gamma, l, gamma_path, "double");
            //TODO peak indices not correct if overlap!=0 until concatenation is
            // not addressed
            save_data(peaks, n_peaks, peaks_path, "int"); // save peak indices
            save_ampd_param(param, param_path);
            save_batch_param(bparam, batch_param_path);
        }
        if(output_lms == 1){
            save_fmtx(lms, lms_path);
        }

        free(data);
        free(sigma);
        free(gamma);
        free(peaks);
        for(int j=0; j<l; j++)
            free(lms->data[j]);
        free(lms->data);
        free(lms);
    }
    free(param);
    free(bparam);
    free(pparam);
    free(conf);
    free(full_data);
    end = clock();
    time_spent = (double)(end - begin) / CLOCKS_PER_SEC; 
    if(verbose > 0)
        printf("runtime = %lf\n",time_spent);
    //fprintf(fp_out, "%d\n", sum_n_peaks);
    fprintf(stdout, "%d\n", sum_n_peaks);
    if(output_peaks == 1)
        fclose(fp_out);
    if(output_rate == 1)
        fclose(fp_out_rate);
    return 0;
}


/**
 * Set data specific hard defined defaults.
 * These can be found in ampd.h. Change accordingly and recompile if needed.
 */
void set_ampd_param(struct ampd_param *p, char *type){

    strcpy(p->datatype,type);
    p->a = DEF_A;
    p->rnd_factor = DEF_RND_FACTOR;
    p->sampling_rate = DEF_SAMPLING_RATE;
    if(strcmp(type, "resp")==0){
        // respiration optimized
        p->sigma_thresh = RESP_SIGMA_THRESHOLD;
        p->peak_thresh = RESP_PEAK_THRESHOLD;

    }
    else if(strcmp(type, "puls")==0){
        // pulsoxy optimized
        p->sigma_thresh = PULS_SIGMA_THRESHOLD;
        p->peak_thresh = PULS_PEAK_THRESHOLD;
    }
    else {
        // default
        p->sigma_thresh = DEF_SIGMA_THRESHOLD;
        p->peak_thresh = DEF_PEAK_THRESHOLD;
    }

}

/**
 * Count occurrences of a character in a file
 *
 */
int count_char(char *path, char cc){
    int count = 0;
    FILE *fp;
    char c;
    fp = fopen(path, "r");
    if(fp == NULL){
        perror("fopen");
        exit(EXIT_FAILURE);
    }
    for (c = getc(fp); c != EOF; c = getc(fp))
        if(c == cc)
            count++;
    fclose(fp);
    return count;
}

/**
 * Extract basename from a path, and omit extension as well.
 * Return 0 on success, -1 on error. Extracted string is put
 * into bname.
 */
void extract_raw_filename(char *path, char *bname, int buf){

    char *tmp;
    memset(bname, 0, buf);
    tmp = basename(path);
    if(strlen(tmp) > buf){
        fprintf(stderr, "extract_raw_filename: insufficent buffer\n");
        exit(EXIT_FAILURE);
    }
    for(int i=0; i<strlen(tmp); i++){
        if(tmp[i] != '.')
            bname[i] = tmp[i];
        else
            break;
    }
    return;
}

/**
 * Function: mkpath
 * ----------------
 *  Recursively create directories. Return -1 on error, 0 on success.
 *  Use mode 0755 for usual read-exec access or 0777 for read-write-exec.
 */
int mkpath(char *file_path, mode_t mode){

    assert(file_path && *file_path);
    for (char* p = strchr(file_path + 1, '/'); p; p = strchr(p + 1, '/')) {
        *p = '\0';
        if (mkdir(file_path, mode) == -1) {
            if (errno != EEXIST) {
                *p = '/';
                return -1;
            }
        }
        *p = '/';
    }
    return 0;
}

/**
 * Save list into file, one value per line.
 * Creates necessary directories.
 * Type represents the type of data: ["float", "double", "int"]
 */
void save_data(void *indata, int n, char *path, char *type){

    FILE *fp;
    float *fdata; int *idata; double *ddata;

    if(mkpath(path, 0777) == -1){
        fprintf(stderr, "cannot make path %s\n",path);
        exit(EXIT_FAILURE);
    }
    fp = fopen(path, "w+");
    if(fp == NULL){
        perror("fopen");
        exit(EXIT_FAILURE);
    }
    // conditional fprint
    if(strcmp(type, "float")==0){
        fdata = (float*)indata;
        for(int i=0; i<n; i++){
            fprintf(fp, "%.3f\n",fdata[i]);
        }
    }
    if(strcmp(type, "double")==0){
        ddata = (double*)indata;
        for(int i=0; i<n; i++){
            fprintf(fp, "%.5lf\n",ddata[i]);
        }
    }
    if(strcmp(type, "int")==0){
        idata = (int*)indata;
        for(int i=0; i<n; i++){
            fprintf(fp, "%d\n",idata[i]);
        }
    }
    fclose(fp);
    return;
}

/**
 * Save matrix to a tab delimited file, without any headers.
 * Creates necessary directories.
 */
void save_fmtx(struct fmtx *mtx, char *path){

    FILE *fp;
    int i, j;
    if(mkpath(path, 0777) == -1){
        fprintf(stderr, "cannot make path %s\n",path);
        exit(EXIT_FAILURE);
    }

    fp = fopen(path, "w");
    if(fp == NULL){
        perror("fopen, exiting...");
        exit(EXIT_FAILURE);
    }
    //printf("rows, cols: %d, %d\n", mtx->rows, mtx->cols);
    for(i=0; i<mtx->rows; i++){
        for(j=0; j<mtx->cols; j++){
            if(j == mtx->cols-1)
                fprintf(fp,"%.3f\n",mtx->data[i][j]);
            else
                fprintf(fp,"%.3f\t",mtx->data[i][j]);
        }
    }
    fclose(fp);
    return;
}
/**
 * Save ampd_param struct contents into file
 */
void save_ampd_param(struct ampd_param *p, char *path){

    FILE *fp;
    fp = fopen(path, "w+");
    if(fp == NULL){
        fprintf(stderr, "cannot open file for writing: %s\n",path);
        exit(1);
    }
    fprintf(fp, "sampling_rate=%lf\n",p->sampling_rate);
    fprintf(fp, "datatype=%s\n", p->datatype);
    fprintf(fp, "a=%lf\n",p->a);
    fprintf(fp, "rnd_factor=%lf\n",p->rnd_factor);
    fprintf(fp, "fit_a=%lf\n",p->fit_a);
    fprintf(fp, "fit_b=%lf\n",p->fit_b);
    fprintf(fp, "fit_r=%lf\n",p->fit_r);
    fprintf(fp, "lambda=%d\n",p->lambda);
    fprintf(fp, "sigma_thresh=%lf\n",p->sigma_thresh);
    fprintf(fp, "peak_thresh=%lf\n",p->peak_thresh);

    fclose(fp);

}
/**
 * Save utility parameters to batch dir
 */
void save_batch_param(struct batch_param *p, char *path){

    FILE *fp;
    fp = fopen(path, "w+");
    if(fp == NULL){
        fprintf(stderr, "can't open file for writing %s\n",path);
        exit(1);
    }
    fprintf(fp,"ind=%d\n",p->ind);
    fprintf(fp,"cycles=%d\n",p->cycles);
    fprintf(fp,"n=%d\n",p->n);
    fprintf(fp,"l=%d\n",p->l);
    fprintf(fp,"batch_length=%lf\n",p->batch_length);
    fprintf(fp,"sampling_rate=%lf\n",p->sampling_rate);
    fprintf(fp,"n_peaks=%d\n",p->n_peaks);
    fprintf(fp,"peaks_per_min=%lf\n",p->peaks_per_min);
    fclose(fp);
}

/**
 * Load a part of the full timeseries data into memory from file.
 * File should only contain one float value on each line.
 */
#define FBUF 32
int fetch_data(char *path, float *data, int n, int ind, int n_zpad){

    FILE *fp;
    int i, count;
    char buf[FBUF];
    char *line = NULL;
    size_t len = 0;
    ssize_t read;
    fp = fopen(path, "r");
    if(fp == NULL){
        perror("fopen");
        exit(EXIT_FAILURE);
    }
    count = 0; i = 0;
    while(fgets(buf, FBUF, fp)){
        if(count < ind){
            count++;
            continue;
        } else if(count > ind +n ){
            break;
        } else {
            sscanf(buf, "%f\n",data + i);
            count++;
            i++;
        }
    }
    fclose(fp);
    return 0;
}
/**
 * Same as fetch_data, but file contents are loaded into memory already
 */
int fetch_data_buff(float *fdata, float *data, int n, int ind, int n_zpad){

    // zpad not used
    int i, count;
    for(i=0;i<n;i++){
        data[i] = fdata[i+ind];
    }
    return 0;

}
/**
 * Load data from file into memory
 */
void load_from_file(char *path, float *fdata, int datalen){

    FILE *fp;
    int i = 0;
    char buf[FBUF];
    fp = fopen(path, "r");
    if(fp == NULL){
        fprintf(stderr, "cannot open file %s\n",path);
        exit(EXIT_FAILURE);
    }
    while(fgets(buf, FBUF, fp) != NULL){
        sscanf(buf, "%f\n",fdata+i);
        i++;
    }
    fclose(fp);
}

void printf_data(float *data, int n){

    for(int i=0; i<n; i++)
        printf("%f\n",data[i]);
    return;
}

void fprintf_data(FILE *fp, float *data, int n){

    for(int i=0; i<n; i++){
        fprintf(fp, "%.5f\n",data[i]);
    }
    return;
}
/*
 * Config handling stuff
 *
 */

void init_preproc_param(struct preproc_param *pparam){
    pparam->preproc = -1; // -1 means unset
    pparam->lpfilt = -1;
    pparam->hpfilt = -1;
}
void set_preproc_param(struct preproc_param *p, char *type){
    if(strcmp(type, "resp") == 0){
        p->preproc = RESP_PREPROC;
        p->hpfilt = RESP_HPFILT;
        p->lpfilt = RESP_LPFILT;
    }
    else if(strcmp(type, "puls")==0){
        p->preproc = PULS_PREPROC;
        p->hpfilt = PULS_HPFILT;
        p->lpfilt = PULS_LPFILT;
    }
    else{
        p->preproc = DEF_PREPROC;
        p->hpfilt = DEF_HPFILT;
        p->lpfilt = DEF_LPFILT;
    }
}

int fork_ampdpreproc(char *path, double lpfilt, double hpfilt){

    return 0;
}
