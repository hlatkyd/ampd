/*
 * filters.h
 *
 * Various filter implementations for smoothing of rat MRI physiological
 * signals, such as respiration, ECG, pulsoxymetry.
 */

void smooth_moving_avg(float *data, int n);

void smooth_savgol(float *data, int n);
