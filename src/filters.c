/*
 * filters.c 
 *
 *
 */


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
 * @param len   Lenght of input data
 * @param m     Halfwidth of window (2*m+1)
 * @param n     Order of fitting
 */
void sgfilt(float *data, int len, int m, int n){

    

}
