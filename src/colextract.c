/*
 * colextract.c
 *
 * Utility program to extract a single column from a delimited data file
 * into another file. Delimiter is found automatically if it is space, comma
 * or tab.
 *
 * Usage from command line:
 * colextract -f [infile] -o [outfile] -n [column index]
 * Optional inputs: -v : print verbose
 *                  -h : print help
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <stdbool.h>
#include <ctype.h>

/**
 * Print general description of input options.
 */
void printf_help(){

    printf(
    "colextract\n"
    "==========\n"
    "Command line utility to extract a column from a delimited data file.\n"
    "The result is copied to a destination file with a single value at\n"
    "each line.\n"
    "\n"
    "Basic usage:\n"
    "colextract -f [inflie] -o [outfile] -n [int]\n"
    "\n"
    "Input options:\n"
    "-h                 print help\n"
    "-f [infile]        path/to/input/file\n"
    "-o [outfile]       path/to/output/file\n"
    "-n [ind]           index of the column to extract\n"
    "-v                 verbose\n"
    );
}
/**
 * Return true if input string represents an integer
 */
bool is_int(char *s){

    int i;
    for (i = 0; i < strlen(s); i++)
    if (isdigit(s[i]) == false)
        return false;

    return true;

}
/**
 * Return the number of occurrences of a character within a string.
 */
int count_char(char *str, char c){

    int count=0;
    for(int i=0; i<strlen(str); i++){
        if(str[i] == c)
            count++;
    }   
    return count;
}
/**
 * Find and return the data delimiter as a string by checking a line.
 */
char get_delim(char *line){

    char d;
    char dlist[] = {'\t',' ',','}; 
    int i, count, count_next;
    count = 0;
    for(i=0; i<sizeof(dlist) / sizeof(dlist[0]); i++){
        count_next = count_char(line, dlist[i]); 
        if(count<count_next){
            count = count_next;
            d = dlist[i];
        }
    }

    return d;
}

#define DELIM "\t" // better not to this
#define MAX_LEN 256
#define BUFS 128
#define VERBOSE_DEFAULT 0
int main(int argc, char **argv){

    int opt;
    int verbose = VERBOSE_DEFAULT;
    char delimchar;
    char delim[2] = {0}; // string delimiter for token
    int n_col = 0;
    char infile[MAX_LEN] = {0};
    char outfile[MAX_LEN] = {0};
    FILE *fp_out;
    // input reading
    FILE *fp_in;
    char *line = NULL;
    size_t len = 0;
    ssize_t read;
    
    char *tok;
    char linebuf[BUFS];
    int checks_done = 0;
    int i, j;
    while((opt = getopt(argc, argv, "ho:f:n:v")) != -1){

        switch(opt){
            case 'n':
                n_col = atoi(optarg);
                break;
            case 'h':
                printf_help();
                return 0;
            case 'o':
                strcpy(outfile, optarg);
                break; 
            case 'f':
                strcpy(infile, optarg);
                break;
            case 'v':
                verbose = 1;
                break;
        }
    }
    if(strcmp(infile,"")==0){
        fprintf(stderr,"No input file given.\n");
        exit(EXIT_FAILURE);
    }
    // print settings
    if(verbose == 1){
        printf("input file: %s\n",infile);
        printf("output file: %s\n",outfile);
        printf("column: %d\n",n_col);
    }
    fp_in = fopen(infile, "r");
    if(fp_in == NULL){
        perror("fopen");
        exit(EXIT_FAILURE);
    }
    fp_out = fopen(outfile, "w");
    if(fp_out == NULL){
        perror("fopen");
        exit(EXIT_FAILURE);
    }
    while((read = getline(&line, &len, fp_in)) != -1){
        if(line[0] == '#' || line[0] == '\n')
            continue; 
        // find delimiters, etc on first line
        if(checks_done == 0){
            delimchar = get_delim(line);
            if(verbose == 1)
                printf("delimiter: '%c'\n", delimchar);
            delim[0] = delimchar;
            delim[1] = '\0';
            checks_done = 1;
        }
        tok = strtok(line, delim);
        if(n_col > 0){
            for(i=0;i<n_col; i++){
                tok = strtok(NULL, delim);
            }
        }
        fprintf(fp_out, "%s\n",tok);


    }
    free(line);
    fclose(fp_in);
    fclose(fp_out);

    return 0;
}

