#!/usr/bin/python3

import matplotlib
matplotlib.use('Qt4Agg')
import matplotlib.pyplot as plt
from pandas.io.parsers import read_csv
import numpy as np
import getopt
import os
import sys

def usage():
    txt="""
Check the results of ampd and saiproc.py by plotting the raw data with the peaks.
Save the figure if needed.

Usage:
    saipcheck.py -i [study_dir]
Optional:
    --save      save png image to study direcotry

The input directory should have the same layout as created by saiproc.py:

[input_dir]/resp.txt
            puls.txt
            ...
            resp.ampd.out/*peaks
            puls.ampd.out/*peaks
            ...
"""
    print(txt)

# data type to plot. names should correspond to the files as in 'usage'
dtype = ["resp", "puls"]

SAVEIMG = False             # save figure in input directory
RAW_OFFSET = 2
PEAKS_OFFSET = 2
NUM = 2                     # currently used number of plots
SAMPLING_RATE = 100         # samples per seconds
RESOLUTION = 1              # only include every nth point to speed up plotting
def main():

    optstr = "hvi:"
    longopt = ["save"]

    if(len(sys.argv) == 1):
        usage()
        sys.exit(0)
    try:
        opts, arg = getopt.getopt(sys.argv[1:],optstr, longopt)
    except getopt.GeoptError as err:
        print(str(err))
        usage()
        sys.exit(2)
    for o, a in opts:
        if o == "-i":
            indir = a;
        elif o == "--save":
            SAVEIMG = True
        else:
            assert False, "unhandled otion"

    fig, axes = plt.subplots(NUM,1,figsize=(14,7))
    for i, d in enumerate(dtype):
        raw = indir +"/"+d + ".txt"
        peaks = indir+"/"+d+".ampd.out/"+d+".peaks"
        #raw_data = np.loadtxt(raw, offset=RAW_OFFSET)
        #peaks_data = np.fromfile(peaks, dtype=np.int, offset=PEAKS_OFFSET)
        raw_data = np.array(read_csv(raw,dtype=float))
        time = np.array([i for i in range(len(raw_data))]) / SAMPLING_RATE
        peaks_x = np.array(read_csv(peaks, dtype=int))
        peaks_x = peaks_x[peaks_x<len(raw_data)]
        peaks_y = raw_data[peaks_x]
        peaks_x = peaks_x / SAMPLING_RATE
        print(peaks_x)
        print("raw_data len="+str(len(raw_data)))
        print("peaks_x max="+str(peaks_x[-1]))
        if d == dtype[0]:
            data_length = len(raw_data)
            max_time = data_length / SAMPLING_RATE

        ax = axes[i]
        ax.plot(time[::RESOLUTION], raw_data[::RESOLUTION],linewidth=0.5)
        ax.scatter(peaks_x, peaks_y, c='r',s=1.5, zorder=9)
        ax.set_aspect('auto')
        ax.set_xlim(left=0,right=max_time)
    plt.tight_layout()
    plt.show()


if __name__ == "__main__":
    main()
