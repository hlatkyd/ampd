/*
 * filters.h
 *
 * Various filter implementations for smoothing of rat MRI physiological
 * signals, such as respiration, ECG, pulsoxymetry.
 */

void moving_avg(float *data, int n);

void sgfilt(float *data, int n);
