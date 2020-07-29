/*
 * rowextract.c
 *
 * Utility to extract samples corresponding to a specific time interval from 
 * a delimited data file which contains columns of uniformly and continously
 * sampled data.
 *
 * Command line usage:
 *
 * rowextract -f [infile] -o [outfile] --start=[n_seconds] --length[n_seconds]
 *
 *
 *
 * Optional inputs:
 *              -h  --help                          print help
 *              -v  --verbose                       verbose
 *              --sampling_rate=[sampling_rate]     in Hz, defaults to 100
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <getopt.h>

// default sampling rate in Hz in case not given as input argument.
// change this if needed, then recompile.
#define DEF_SAMPLING_RATE 100

// copy comment or header lines from the beginning of input file
#define COPY_HEADER 1

#define MAX_PATH_LEN 1024
static struct option long_options[]=
{
    {"infile",required_argument, NULL, 'f'},
    {"outfile",required_argument, NULL, 'o'},
    {"start",required_argument, NULL, 's'},
    {"length",required_argument,NULL,'l'},
    {"sampling-rate",required_argument,NULL,'r'},
    {"verbose",no_argument,NULL,'v'},
    {"help",no_argument, NULL,'h'},
    {0,0,0,0},
};
static char optstring[] = "hvf:o:s:l:r:";
void printf_help(){

    printf(
     "\n"
     "rowextract\n"
     "==========\n"
     "Utility to extract samples corresponding to a specific time interval from\n "
     "a delimited data file which contains columns of uniformly and continously\n"
     "sampled data.\n"

     "Command line usage:\n"
     "\n"
     "rowextract -f [infile] -o [outfile] --start=[n_seconds] --length[n_seconds]\n"
     "\n"
     "Required inputs:\n"
     "-f    --infile                    input file path\n"
     "-o    --outfile                   output file path\n"
     "-s    --start                     start time of extracted interval in seconds\n"
     "-l    --length                    length of extracted interval in seconds\n"
     "Optional inputs:\n"
     "-h    --help                      print help\n"
     "-v    --verbose                   verbose\n"
     "--sampling-rate=[sampling_rate]   in Hz, defaults to 100\n"
     
    );
}

int main(int argc, char **argv){

    // input arguments
    int verbose = 0;
    char infile[MAX_PATH_LEN] = {0};
    char outfile[MAX_PATH_LEN] = {0};
    FILE *fp_in, *fp_out;
    double sampling_rate = DEF_SAMPLING_RATE;
    double start_time = -1;
    double length = -1;
    // for reading lines
    char *line = NULL;
    size_t len = 0;
    ssize_t read;
    // cope input file header to output
    int copy_header = COPY_HEADER;
    int count = 0;
    int minind, maxind; // indices of first and last line
    int opt;
    while((opt = getopt_long(argc, argv, optstring, long_options, NULL)) != -1){
        switch(opt){
            case 'v':
                verbose = 1;
                break;
            case 'h':
                printf_help();
                return 0;
            case 's':
                start_time = atof(optarg);
                break;
            case 'l':
                length = atof(optarg);
                break;
            case 'r':
                sampling_rate = atof(optarg);
                break;
            case 'f':
                strcpy(infile, optarg);
                break;
            case 'o':
                strcpy(outfile,optarg);
                break;
        }
    }
    if(argc == 1){
        printf_help();
        return 0;
    }
    // checking if inputs are sufficent
    if(length == -1){
        fprintf(stderr,"Please specify interval length: -l [len] --length=[len]\n");
        exit(EXIT_FAILURE);
    }
    if(start_time == -1){
        fprintf(stderr,"Please specify interval start: -s [seconds] --start=[seconds]\n");
        exit(EXIT_FAILURE);
    }
    if(strcmp(infile,"") == 0){
        fprintf(stderr,"Please specify input file.\n");
        exit(EXIT_FAILURE);
    }
    if(strcmp(outfile,"") == 0){
        fprintf(stderr,"Please specify output file.\n");
        exit(EXIT_FAILURE);
    }

    //
    minind = (int)(start_time * sampling_rate);
    maxind = (int)((start_time + length) * sampling_rate);
    fp_in = fopen(infile,"r");
    fp_out = fopen(outfile, "w+");
    if(fp_in == NULL){
        fprintf(stderr, "Cannot open file on path: %s\n",infile);
        exit(EXIT_FAILURE);
    }
    if(fp_out == NULL){
        fprintf(stderr, "Cannot open file on path: %s\n",outfile);
        exit(EXIT_FAILURE);
    }
    if(verbose == 1){
        printf("time interval: %lf sec, sampling rate: %lf\n",length,sampling_rate);
        printf("Copying rows %d - %d to file %s\n",minind, maxind,outfile);
    }
    while((read = getline(&line,  &len, fp_in)) != -1){
        if(line[0] == '#' || line[0] == '\n'){
            if(copy_header == 1 && line[0] == '#'){
                fprintf(fp_out,"%s",line);
            }
            continue;
        }
        else if(count < minind){
            count++;
            continue;
        }
        else if(maxind < count){
            break;
        }
        // copy line to output file
        else{
            fprintf(fp_out,"%s",line);
            count++;
        }

    }
    fclose(fp_in);
    fclose(fp_out);

    return 0;
}
