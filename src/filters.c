/*
 * filters.c 
 *
 *
 */

#include "filters.h"

/**
 * Simple low pass filter in time domain as per wiki page
 *
 * @param data          input timeseries
 * @param n             length of data
 * @param sample_rate   sampling rate in seconds
 * @param cutoff_freq   cutoff frequency in Hz
 *
 */
void tdlpfilt(float *data, int n, double sample_rate, double cutoff_freq){

    int i;
    float rc = 1.0 / (cutoff_freq * 2*3.14);
    float dt = 1.0 / sample_rate;
    float alpha = dt / (rc + dt);
    float *buf = malloc(sizeof(float) * n);
    
    for(i=0; i<n; i++){
        buf[i] = data[i];
    }
    for(i=1; i<n; i++){
        data[i] = data[i-1] + alpha * (buf[i] - data[i-1]);
    }
    free(buf);
}

/**
 * Simple high pass filter in time domain as per wiki page
 *
 * @param data          input timeseries
 * @param n             length of data
 * @param sample_rate   sampling rate in seconds
 * @param cutoff_freq   cutoff frequency in Hz
 *
 */
void tdhpfilt(float *data, int n, double sample_rate, double cutoff_freq){

    int i;
    float rc = 1.0 / (cutoff_freq * 2*3.14);
    float dt = 1.0 / sample_rate;
    float alpha = rc / (rc + dt);
    float *buf = malloc(sizeof(float) * n);
    
    for(i=0; i<n; i++){
        buf[i] = data[i];
    }
    for(i=1; i<n; i++){
        data[i] = alpha * (data[i-1] + buf[i] - buf[i-1]);
    }
    free(buf);
}
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
