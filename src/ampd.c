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
int output_meta = DEF_OUTPUT_META;
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
    "\nampd - Automatic Multi-scale Peak detection - version %d.%d\n"
    "=========================================================================\n"
    ,VERSION_MAJ,VERSION_MIN);
    printf(
    "Peak detection algorithm for quasiperiodic data, as specified in paper:\n"
    "An Efficient Algorithm for Automatic Peak Detection in Noisy Periodic and\n"
    "Quasi-Periodic Signals. doi:10.3390/a5040588\n"
    "This program takes a file input which contains a single timeseries of\n"
    "quasi-periodic data. The main output is a direcory with files containing\n"
    "the peak indices, peak rates and metadata. An aux output directory can be\n"
    "created as well for inspection of all intermediary data processing.\n"
    "The final peak count is sent to stdout.\n"
    "\n"
    "The input file should only contain a single series of data, with 1 value\n"
    "each line.\n"

    "\n"
    "Usage from command line:\n"
    "ampd -f [infile] -o [out_dir] -r [sampling_rate] -l [batch_length]\n"
    "\n"
    "Optional arguments:\n"
    "-o --outdir:           output dir, defaults to [cwd]/ampd.out\n"
    "-v --verbose:          verbose\n"
    "-h --help:             print help\n"
    "-r --samplig-rate:     sampling rate input data in Hz, default is 100\n"
    "-l --batch-length:     data window length in seconds, default is 60 sec\n"
    "--preproc:             call ampdpreproc on data first\n"
    "--lpfilt:              apply lowpass filter in Hz\n"
    "--hpfilt:              apply highpass filter in Hz\n"
    "-a --auxdir:           aux data root dir, default is [cwd]/ampd.aux\n"
    "--output-all:          output aux data, local maxima scalogram not included\n"
    "--output-lms:          output local maxima scalogram (high disk space usage)\n"
    "--output-rate:         output peak-per-min\n"
    "--output-peaks         output peak indices\n"
    "\n"
        );
}
/**
 * Main handles the command line input parsing and output file management.
 */
int main(int argc, char **argv){

    int i, j, opt;
    struct ampd_config *conf;
    char datatype[32] = {0};
    // timing
    clock_t begin, end;
    double time_spent;

    char infile[MAX_PATH_LEN] = {0};
    char infile_basename[MAX_PATH_LEN]; // input, without dir and extension
    char outdir[MAX_PATH_LEN] = {0};
    char outfile_peaks[MAX_PATH_LEN] = {0}; // main output file with indices of peaks
    char outfile_rate[MAX_PATH_LEN] = {0};  // rate per min for each batch
    char outfile_meta[MAX_PATH_LEN] = {0}; // metadata
    FILE *fp_out;        // main output file containing the peak indices
    FILE *fp_out_rate;
    FILE *fp_out_meta;
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
    double batch_length = -1;
    int cycles;         // number of data batches
    int sum_n_peaks;    // summed peak number from all batches
    int *sum_peaks;      // concatenated vector from peak indices
    double peaks_per_min;
    struct batch_param *bparam; // only for outputting batch utility parameters

    /* 
     * filtering
     */
    struct preproc_param *pparam;
    int preproc = -1;
    double hpfilt = -1;
    double lpfilt = -1;
    /*
     * ampd routine
     */
    struct ampd_param *param;
    double sampling_rate = -1;
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
    struct meta_param *mparam;

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
                break;
            case 'o':
                strcpy(outdir, optarg);
                break;
            case 'f':
                strcpy(infile, optarg);
                break;
            case 'r':
                sampling_rate = atof(optarg);
                break;
            case 'a':
                strcpy(aux_dir, optarg);
                break;
            case 'l':
                batch_length = atof(optarg);
                break;
            case ARG_OUTPUT_ALL:
                output_all = 1;
                break;
            case ARG_OUTPUT_LMS:
                output_lms = 1;
                break;
            case ARG_OUTPUT_PEAKS:
                output_peaks = 1;
                break;
            case ARG_OUTPUT_RATE:
                output_rate = 1;
                break;
            case ARG_PREPROC:
                // do filtering and smoothing
                preproc = 1;
                break;
            case ARG_HPFILT:
                hpfilt = atof(optarg);
                break;
            case ARG_LPFILT:
                lpfilt = atof(optarg);
                break;
            case ARG_OUTPUT_IMG:
                output_img = 1;
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
    if(strcmp(infile,"")==0){
        fprintf(stderr, "No input file specified.\n");
        //free_conf_malloc_onerr();
        exit(EXIT_FAILURE);
    }
    extract_raw_filename(infile, infile_basename, sizeof(infile_basename));
    if(strcmp(datatype,"")==0){
        strcpy(datatype,"def");
        //strcpy(conf->datatype,datatype);
    }
    set_ampd_param(param, datatype);
    set_preproc_param(pparam, datatype);
    if(preproc != -1)
        pparam->preproc = preproc;
    if(lpfilt != -1)
        pparam->lpfilt = lpfilt;
    if(hpfilt != -1)
        pparam->hpfilt = hpfilt;
    if(strcmp(outdir,"")==0){
        snprintf(outdir,sizeof(outdir),"%s/%s",cwd,outdir_def);
    }

    if(strcmp(aux_dir,"")==0){
        snprintf(aux_dir, sizeof(aux_dir), "%s/%s",cwd, aux_dir_def);
    }
    //TODO fix
    if(sampling_rate == -1)
        sampling_rate = DEF_SAMPLING_RATE;
    if(batch_length == -1)
        batch_length = DEF_BATCH_LENGTH;

    param->sampling_rate = sampling_rate;
    // setting outptu files
    snprintf(outfile_peaks,sizeof(outfile_peaks),"%s/%s.peaks",outdir,infile_basename);
    snprintf(outfile_rate,sizeof(outfile_rate),"%s/%s.rate",outdir,infile_basename);
    snprintf(outfile_meta,sizeof(outfile_meta),"%s/%s.meta",outdir,infile_basename);
    // setting available param
    // set available config
    // setting remaining variables for processing
    sum_n_peaks = 0;
    datalen = count_char(infile, '\n');
    data_buf = (int)(batch_length * param->sampling_rate);
    cycles = (int) (ceil(datalen / (double) (data_buf)));

    /* fill batch param */
    bparam->cycles = cycles;
    bparam->batch_length = batch_length;
    bparam->sampling_rate = sampling_rate;

    // preload data and apply filters
    full_data = malloc(sizeof(float) * datalen);
    load_from_file(infile, full_data, datalen);
    // make path
    // opening main output files
    if(output_peaks == 1){
        mkpath(outfile_peaks, 0777);
        fp_out = fopen(outfile_peaks, "w+");
        if(fp_out == NULL){
            fprintf(stderr, "cannot open file for writing %s\n",outfile_peaks);
            exit(EXIT_FAILURE);
        }
    }
    if(output_rate == 1){
        mkpath(outfile_rate, 0777);
        fp_out_rate = fopen(outfile_rate, "w+");
        if(fp_out_rate == NULL){
            fprintf(stderr,"cannot open file %s\n",outfile_rate);
            exit(EXIT_FAILURE);
        }
        fprintf(fp_out_rate,"# batch_length=%lf\n",batch_length);
        fprintf(fp_out_rate,"# sampling_rate=%lf\n",param->sampling_rate);

    }
    if(verbose > 0){
        printf("ampd input\n-----------------\n");
        printf("verbose: %d\n",verbose);
        printf("infile: %s\n", infile);
        printf("outdir: %s\n", outdir);
        printf("output-all: %d\n",output_all);
        printf("aux_dir: %s\n", aux_dir);
        printf("datatype: %s\n",param->datatype);
        printf("preproc=%d\n",pparam->preproc);
        printf("lowpassfilt=%lf\n",pparam->lpfilt);
        printf("highpassfilt=%lf\n",pparam->hpfilt);
        printf("sampling_rate: %.5lf\n",param->sampling_rate);
        printf("batch_length: %lf\n",batch_length);
        printf("datalen: %d\n", datalen);
        printf("data_buf: %d\n",data_buf);
        printf("cycles: %d\n", cycles);
        printf("output-lms: %d\n",output_lms);
        printf("output-rate: %d\n",output_rate);

    }

    // malloc
    n = (int)data_buf;
    l = (int)ceil(n/2)-1;
    bparam->n = n;
    bparam->l = l;
    lms = malloc_fmtx(l, n);
    data = malloc(sizeof(float)*n);
    gamma = malloc(sizeof(double)*l);
    sigma = malloc(sizeof(double)*n);
    peaks = malloc(sizeof(int)*n);

    /*
     * Processing
     *
     */
    for( i=0; i<cycles; i++){
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

        ind = i * data_buf;
        bparam->ind = ind;

        if(verbose > 1){
            printf("\nfetchig data:\n");
            printf("infile=%s\n",infile);
            printf("data=%p\n",data);
            printf("n=%d, ind=%d\n",n,ind);
        }
        // load data batch
        fetch_data_buff(full_data, datalen, data, n, ind, n_zpad);
        if(output_all == 1)
            save_data(data, n, raw_path,"float"); // save raw data

        // preproc
        linear_fit(data, n, param);
        linear_detrend(data, n, param);
        if(output_all == 1)
            save_data(data, n, preproc_path,"float"); // save detrend data
        if(pparam->preproc == 1){
            if(pparam->hpfilt > 0){
                tdhpfilt(data, n, param->sampling_rate, pparam->hpfilt);
            }
            if(pparam->lpfilt > 0){
                tdlpfilt(data, n, param->sampling_rate, pparam->lpfilt);
            }
        }

        // main ampd routine, contains malloc
        n_peaks = ampdcpu(data, n, param, lms, gamma, sigma, peaks);

        //catch_false_pks(peaks, &n_peaks, ts, thresh);
        sum_n_peaks += (int)(n_peaks );
        if(verbose > 0){
            printf("batch=%d/%d, n=%d, sum=%d, "
                    "mean_dst=%.3lf s, stdev_dst=%.3lf s\n",
                    i,cycles, n_peaks,sum_n_peaks,
                    param->mean_pk_dist, param->stdev_pk_dist);
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
            save_data(peaks, n_peaks, peaks_path, "int"); // save peak indices
            save_ampd_param(param, param_path);
            save_batch_param(bparam, batch_param_path);
        }
        if(output_lms == 1){
            save_fmtx(lms, lms_path);
        }
    }
    if(output_peaks == 1)
        fclose(fp_out);
    if(output_rate == 1)
        fclose(fp_out_rate);
    // save some metadata to file
    if(output_meta == 1){
        mparam = malloc(sizeof(struct meta_param));
        memset(mparam, 0, sizeof(mparam));
        strcpy(mparam->infile, infile);
        strcpy(mparam->basename, infile_basename);
        strcpy(mparam->datatype, datatype);
        mparam->sampling_rate = bparam->sampling_rate;
        mparam->batch_length = bparam->batch_length;
        mparam->total_peaks = sum_n_peaks;
        mparam->total_batches = cycles;
        save_meta(mparam, pparam,  outfile_meta);
        free(mparam);
    }
    // free ampd data arrays
    free(data);
    free(sigma);
    free(gamma);
    free(peaks);
    for(j=0; j<l; j++)
        free(lms->data[j]);
    free(lms->data);
    free(lms);

    // free parameters and stuff
    free(param);
    free(bparam);
    free(pparam);
    free(conf);
    free(full_data);
    // finalize
    end = clock();
    time_spent = (double)(end - begin) / CLOCKS_PER_SEC; 
    if(verbose > 0)
        printf("runtime = %lf sec\n",time_spent);
    return sum_n_peaks;
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
    int i;
    memset(bname, 0, buf);
    tmp = basename(path);
    if(strlen(tmp) > buf){
        fprintf(stderr, "extract_raw_filename: insufficent buffer\n");
        exit(EXIT_FAILURE);
    }
    for(i=0; i<strlen(tmp); i++){
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
    int i;
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
        for(i=0; i<n; i++){
            fprintf(fp, "%.3f\n",fdata[i]);
        }
    }
    if(strcmp(type, "double")==0){
        ddata = (double*)indata;
        for(i=0; i<n; i++){
            fprintf(fp, "%.5lf\n",ddata[i]);
        }
    }
    if(strcmp(type, "int")==0){
        idata = (int*)indata;
        for(i=0; i<n; i++){
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
    fprintf(fp, "mean_pk_dist=%.3lf\n",p->mean_pk_dist);
    fprintf(fp, "stdev_pk_dist=%.3lf\n",p->stdev_pk_dist);

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
    //fprintf(fp, "mean_pk_dist=%.3lf\n",p->mean_pk_dist);
    //fprintf(fp, "stdev_pk_dist=%.3lf\n",p->stdev_pk_dist);
    fclose(fp);
}
void save_meta(struct meta_param *p, struct preproc_param *pp, char *path){

    FILE *fp;
    mkpath(path, 0777);
    fp = fopen(path, "w+");
    if(fp == NULL){
        fprintf(stderr, "cannot open file %s\n",path);
        exit(EXIT_FAILURE);
    }
    fprintf(fp,"infile=%s\n",p->infile);
    fprintf(fp,"basename=%s\n",p->basename);
    fprintf(fp,"datatype=%s\n",p->datatype);
    fprintf(fp,"preproc=%d\n",pp->preproc);
    fprintf(fp,"hpfilt=%lf\n",pp->hpfilt);
    fprintf(fp,"lpfilt=%lf\n",pp->lpfilt);
    fprintf(fp,"sampling_rate=%lf\n",p->sampling_rate);
    fprintf(fp,"batch_length=%lf\n",p->batch_length);
    fprintf(fp,"total_batches=%d\n",p->total_batches);
    fprintf(fp,"total_peaks=%d\n",p->total_peaks);
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
 * Same as fetch_data, but file contents are loaded into memory already.
 * n number of values are loaded starting frim index 'ind'. If it would exceed
 * the length of original, the last n values are loaded regarless of ind.
 */
#define IND_READJUST 1
int fetch_data_buff(float *fdata,int len, float *data, int n, int ind, int n_zpad){

    // zpad not used
    int i, count;
    if(IND_READJUST == 1){
        if(ind + n <len){
            for(i=0;i<n;i++){
                data[i] = fdata[i+ind];
            }
        }
        else{
            for(i=0;i<n;i++){
                data[i] = fdata[len-n+i];
            }
        }
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

    int i;
    for(i=0; i<n; i++)
        printf("%f\n",data[i]);
    return;
}

void fprintf_data(FILE *fp, float *data, int n){

    int i;
    for(i=0; i<n; i++){
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
/**
 * Set preprocessing parameters based on datatype to defaults.
 *
 * Modify these parameters later on.
 */
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
