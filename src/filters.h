/*
 * filters.h
 *
 * Various filter implementations for smoothing of rat MRI physiological
 * signals, such as respiration, ECG, pulsoxymetry.
 */

#include <stdlib.h>
#include <stdio.h>
#include <math.h>

void movingavg(float *data, int n, int w);
void hpfilt(float *data, int n, double sample_rate, double cutoff_freq);

void sgfilt(float *data, int n, int w, int p);
