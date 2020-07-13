/*
 * filters.c 
 *
 *
 */

#include "filters.h"

/**
 * Apply moving average smoothing to floating point uniformly sampled data.
 * End points are not truncated, but weighted. The resultant data is put in
 * same memory space as the input.
 *
 * @param data      input data
 * @param n       length of input data
 * @param w         half of averaging window (2*w+1)
 */
void movingavg(float *data, int n, int w){

    float *buf;
    int i, j;
    /* from data x0 x1 x2 x3 x4 x5 ...
     * make buffer of n+2*w length
     * xw x(w-1) ... x1 x0 x1 x2 x3 ...
     *
     */
    printf("w=%d\n",w);
    buf = malloc(sizeof(float)*(n+2*w));
    for(i=0; i<n; i++)
        buf[i+w] = data[i];
    // fill  end points for buffer
    for(i=0; i<w; i++){
        buf[i] = data[w-i]; // start
        buf[n+i] = data[n-i];  // end
    }
    for(i=0; i<n; i++){
        data[i] = 0.0;
        for(j=0; j<w*2+1; j++){
            data[i] += buf[i+j]/(2*w+1);
        }
        printf("%lf\n",data[i]);
    }
    for(i=0; i<n; i++){
        //printf("%lf\n",buf[i]);
    }
    free(buf);
}

/**
 * Apply Savicky-Golay filter to uniformly sampled input data of length n,
 * by use of Gram-polynomials. End points are treated without truncation. 
 * 
 * Ref: General Least-Squares Smoothing and Differentiation by the Convolution
 * (Savitzky-Golay) Method
 *
 * https://pubs.acs.org/doi/pdf/10.1021/ac00205a007
 * The result is saved in place of the input data.
 *
 * @param data  Uniformly sampled data
 * @param n   Lenght of input data
 * @param w     Halfwidth of window (2*m+1)
 * @param p     Order of fitting
 */
void sgfilt(float *data, int n, int w, int p){

    

}
