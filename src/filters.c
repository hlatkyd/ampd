/*
 * filters.c 
 *
 *
 */


/**
 * Apply moving average smoothing to floating point uniformly sampled data.
 * End points are not truncated, but weighed. The resultant data is put in
 * same memory space as the input.
 *
 * @param data      input data
 * @param len       length of input data
 * @param w         half of averaging window (2*w+1)
 */
void movingavg(float *data, int len, int w){

    float *buf;
    float va;
    int i, j;
    buf = malloc(sizeof(float) * (2*w+1));
    for(i=0; i<len; i++){
        // beginning points
        if(i<w+1){
            for(j=-m; j< 2*m; j++){

            } 
        }
    }
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
 * @param len   Lenght of input data
 * @param m     Halfwidth of window (2*m+1)
 * @param n     Order of fitting
 */
void sgfilt(float *data, int len, int m, int n){

    

}
