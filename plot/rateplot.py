#!/usr/bin/python3

"""
Script to plot physiological data made with AMPD

"""
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import FormatStrFormatter
import sys
import os
import glob
from pathlib import Path    # Python > 3.5 
import csv

def usage():

    txt="""



    """
    print(txt)


# USER SETUP
#------------------------------------------------------------------------------
# Directory to contain processed physiological data for each study.
STUDY_ROOT_DIR="/home/david/work/proc"
# Experiment log file
STUDY_LOG = "/home/david/work/jk2020aug.csv"
# Header file for expriment log
STUDY_LOG_HEADER = "/home/david/work/csv_headers"
# File to list included studies i the plotting and analysis
STUDY_ID_FILE = "physplot_scop_id"

# plots
FIG_SIZE=(12,8)

# Study directories are searched as $IS_PREFIX[YYYYMMDDNN]$ID_SUFFIX
ID_PREFIX = "s_"    # prefix for study directory
ID_SUFFIX = ""      # suffix for study directory
# AMPD rate processing
HEADER_LENGTH = 2   # maximum number of header lines in .rate files
# Log file processing
TIMESTEP = 60       # pick points

#------------------------------------------------------------------------------

def main():

    pass


def read_rate_files(study_list, dtype=None):
    """Return  data array, and metadata as batch length, sampling rates"""
    
    if dtype == "resp":
        glob_pattern = "*resp.rate"
    elif dtype == "puls":
        glob_pattern = "*puls.rate"

    indir = full_path(STUDY_ROOT_DIR)
    #---------- gather files to plot-------------------------------------------
    pattern = ID_PREFIX + "*"+ID_SUFFIX
    study_dir_list = []
    data_files = []

    for studydir in sorted(glob.glob(indir+"/"+pattern)):

        if os.path.basename(studydir) in study_list:
            study_dir_list.append(studydir)
            # search files
            for path in Path(studydir).rglob(glob_pattern):
                if DTYPE[i] in str(path.name):
                    data_files[i].append(path)

    # --------------------read files-------------------------------------------
    data_arr = [[] for i in range(NUM)]
    batch_lengths = []
    sampling_rates = []
    for ind, data_file in enumerate(data_files):
        for f in data_file:
            with open(f,'r') as openf:
                data = openf.read().split('\n')
                for item in data[:HEADER_LENGTH]:
                    if "batch_length" in item:
                        batch_length = int(float(item.split('=')[-1]))
                        batch_lengths.append(batch_length)
                    if "sampling_rate" in item:
                        val = int(float(item.split('=')[-1]))
                        sampling_rates.append(val)
                data = np.array([int(x) for x in data[:-1] if not x[0] == "#"])
                data_arr[ind].append(data)

    # check for samebatch lengts and sampling rates
    bl = batch_lengths[0]
    sr = sampling_rates[0]
    if any(t != bl for t in batch_lengths):
        print("WARNING: Batch lengths are not the same.")
    if any(t != sr for t in sampling_rates):
        print("WARNING: Sampling rates are not the same.")

    return data_arr, batch_lengths, sampling_rates


if __name__ == "__main__":
    main()
