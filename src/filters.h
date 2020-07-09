/*
 * filters.h
 *
 * Various filter implementations for smoothing of rat MRI physiological
 * signals, such as respiration, ECG, pulsoxymetry.
 */

void movingavg(float *data, int n, int w);

void sgfilt(float *data, int n);
