#!/usr/bin/python3

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import FormatStrFormatter
import sys
import os
import glob
from pathlib import Path    # Python > 3.5 
import csv
"""
Script to plot physiological data



"""
# USER SETUP
#------------------------------------------------------------------------------
# Directory to contain processed physiological data for each study.
STUDY_ROOT_DIR="/home/david/work/saidata_proc"
# Experiment log file
STUDY_LOG = "/home/david/work/jk2020aug.csv"
# Header file for expriment log
STUDY_LOG_HEADER = "/home/david/work/csv_headers"
# File to list included studies i the plotting and analysis
STUDY_ID_FILE = "physplot_study_id"
#STUDY_ID_FILE = "physplot_scop_id"
# plots
PLOT_INDIVIDUAL_RATES = True
PLOT_AVG_RATES = True


#------------------------------------------------------------------------------

# DEV SETUP
#------------------------------------------------------------------------------
# Study directories are searched as $IS_PREFIX[YYYYMMDDNN]$ID_SUFFIX
ID_PREFIX = "s_"    # prefix for study directory
ID_SUFFIX = ""      # suffix for study directory

# AMPD rate processing
# --------------------
DTYPE = ["resp","puls"]
DTYPE_THRESH = [(30,150),(100,500)]
NUM = len(DTYPE)    # number of data type, plots
HEADER_LENGTH = 2   # maximum number of header lines in .rate files

# Log file processing
# -------------------
TIMESTEP = 60       # pick points

#------------------------------------------------------------------------------



def usage():

    txt="""


    """
    print(txt)

def main():

    print("study_id file: "+str(STUDY_ID_FILE))
    print("study log file: "+str(STUDY_LOG))

    study_list = read_study_id_file(STUDY_ID_FILE)

    data, batch_length, sampling_rate = read_rate_files(study_list)

    sdict = fetch_study_log(study_list, STUDY_LOG, STUDY_LOG_HEADER)

    #rateplot(study_list, data, batch_length[0], sampling_rate[0])

    logplot(sdict)

    plt.show()
    return 0

# MAIN FUNCTIONS

def logplot(sdict):
    """ Plot data from experiment log file
    
    Missing instances are plotted with dotted line
    
    """
    YLIM_ISO=(0,1)
    YLIM_RESP=(20,120)
    YLIM_PULS=(100,500)
    FILL_GAPS=False
    """
    print(sdict[0].keys())
    print(type(sdict))
    print(type(sdict[0]))
    print(sdict[0]["isoflurane"])
    print(len(sdict))
    """
    # check consistency

    # find longest study
    # generate x axis time points
    maxtime_m = int(int(max([x["time"][-1] for x in sdict]))/TIMESTEP)
    data_x = np.arange(0,maxtime_m,1)

    # create y data
    for sd in sdict:
        time = [int(int(x)/TIMESTEP) if x is not '' else False for x in sd["time"]]
        sd["time_m"] = np.array(time)

        # generate lines
        sd["iso_m"] = np.array([float(x) if x is not '' else np.nan for x in sd["isoflurane"]])
        sd["puls_m"] = np.array([float(x) if x is not '' else np.nan for x in sd["bpm"]])
        sd["resp_m"] = np.array([float(x) if x is not '' else np.nan for x in sd["resp"]])

    # plot
    num = len(sdict)
    y_subplots = 4
    x_subplots = int(num / y_subplots)
    if num % y_subplots != 0:
        x_subplots += 1
    fig, axes = plt.subplots(x_subplots,y_subplots,figsize=(15,7))
    fig.subplots_adjust(wspace=0.5)
    for i, ax in enumerate(axes.flatten()):
        try:
            x = sdict[i]["time_m"]
            #y_iso = interpolate_gaps(sdict[i]["iso_m"], limit=2)
            y_iso = sdict[i]["iso_m"]
            y_puls = sdict[i]["puls_m"]
            y_resp = sdict[i]["resp_m"]
            if FILL_GAPS:
                fill_gaps(x, dtype="time")
                fill_gaps(y_iso)
                fill_gaps(y_puls)
                fill_gaps(y_resp)

                pass
            ax_resp = ax.twinx()
            ax_puls = ax.twinx()
            ax.plot(x, y_iso, c='r')
            ax_puls.plot(x, y_puls,c='g')
            ax_resp.plot(x, y_resp,c='b')
            ax.set_ylim(YLIM_ISO)
            ax_resp.set_ylim(YLIM_RESP)
            ax_puls.set_ylim(YLIM_PULS)
            ax.tick_params(axis='y',colors='r')
            ax_resp.tick_params(axis='y',direction='out',pad=25, colors='b')
            ax_puls.tick_params(axis='y', colors='g')
        except Exception as e:
            print(e)

    pass

def rateplot(study_list, data_arr, batch_lengths, sampling_rates):
    """ """
    # plot individual lines

    #study_list = [os.path.basename(x) for x in study_dir_list]
    
    fig, axes = plt.subplots(NUM,1, figsize=(15,7))
    data_xmax = 0
    # plot average
    if PLOT_AVG_RATES:
        data_arr_2 = [0 for i in range(len(DTYPE))]
        for i in range(len(DTYPE)):
            minlen = min([len(data) for data in data_arr[i]])
            maxlen = max([len(data) for data in data_arr[i]])
            data_arr_2[i] = [data[:minlen] for data in data_arr[i]]

    # make avg
    arr = np.array(data_arr_2)
    avg = arr.mean(axis=1)
    std = np.std(arr, axis=1)
    for num, ax in enumerate(axes):
        #for n, data in enumerate(data_list[num]):
        data = avg[num]
        data_stddev = std[num]
        ax.plot(data, color="black", label="avg", linewidth=2,zorder=100)
        ci1 = data + 2 * data_stddev
        ci2 = data - 2 * data_stddev
        x = np.arange(0, minlen)
        ax.fill_between(x, ci1, ci2, color="blue", alpha=0.1)

    for num, ax in enumerate(axes):
        for n, data in enumerate(data_arr[num]):
            ax.plot(data, label=study_list[n], alpha=0.8)
            if len(data) > data_xmax:
                data_xmax = len(data)
        maj_tick_freq = batch_lengths
        min_tick_freq = maj_tick_freq / 6
        ax.set_xlim(left=0,right=data_xmax)
        major_ticks = np.arange(0,data_xmax, maj_tick_freq)
        minor_ticks = np.arange(0,data_xmax, min_tick_freq)
        ax.set_xticks(major_ticks)
        ax.set_xticks(minor_ticks,minor=True)
        ax.tick_params(axis='both',which='major',labelsize=10)
        ax.xaxis.set_major_formatter(FormatStrFormatter('%d min'))
        ax.tick_params(axis='both',which='minor',labelsize=0)
        #ax.xaxis.grid(True, which='minor')
        pos1 = ax.get_position()
        pos2 = [0.05, pos1.y0, pos1.width*0.9, pos1.height]
        ax.set_position(pos2)
        ax.set_title(label=str(DTYPE[num]))
        if num == 0:
            # generate text from excluded studies
            props = dict(boxstyle='square', facecolor='white', alpha=0.5)
            ax.legend(loc='upper left',bbox_to_anchor=(1.,1))



def read_rate_files(study_list):
    """Return  data array, and metadata as batch length, sampling rates"""

    indir = full_path(STUDY_ROOT_DIR)
    #---------- gather files to plot-------------------------------------------
    pattern = ID_PREFIX + "*"+ID_SUFFIX
    study_dir_list = []
    data_files = [[] for i in range( NUM )]

    for studydir in sorted(glob.glob(indir+"/"+pattern)):

        if os.path.basename(studydir) in study_list:
            study_dir_list.append(studydir)
            # search files
            for path in Path(studydir).rglob('*.rate'):
                for i in range(NUM):
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




# UTIL FUNCTIONS

def read_study_id_file(path):
    """ Return studyid list from study id file"""
    path = full_path(path)
    if not os.path.isfile(path):
        print("Wrong path for study_id file: '"+str(path)+"'")
        sys.exit(0)
    study_list = []
    with open(path) as f:
        lines = f.readlines()
        for line in lines:
            if line[0] in ("#","\t"," ","\n"):
                continue
            else:
                study_list.append(line.split("\n")[0])
    # check format
    for s in study_list:
        if s.startswith(ID_PREFIX) and s.endswith(ID_SUFFIX):
            continue
        else:
            print("ERROR: wrong study string format in file "+str(path))
            sys.exit(0)
    return sorted(study_list)

def full_path(path):
    """ check path syntax and return full path"""
    if path[0] == "~":
        return os.path.expanduser("~")+"/"+path[1:]
    else:
        return os.path.abspath(path)

def fetch_study_log(study_list, jk_path, header_path):
    """
    Return a list a dictionaries, one dictionary per study.
    One dictionary contains lists of various parameters, one element 
    corresponding to a sequence.
    Keys are found in jk header file 'csv_headers' under 'full_log'

    currently the keys are:

    ['studyid', 'pslabel', 'comment', 'scantime', 'tr', 'te', 'images', 'nt',
    'ss', 'seqid', 'fwhm', 'type', 'ratid', 'weight', 'measurement',
    'anesthesia', 'treatment', 'time', 'resp', 'bpm', 'isoflurane',
    'resp_var', 'bpm_var', 'op_expcomment', 'a_expcomment',
    'Prescan_FatOffset', 'H1offset', 'pwr90']

    """

    study_dict_list = []
    study_row_list = [[] for i in range(len(study_list))]
    header_list = read_csv_headers(header_path)

    # read csv file and sort relevant rows into list of lists
    with open(jk_path) as csv_file:
        csv_reader = csv.reader(csv_file,delimiter=',')
        line_count = 0
        curstudy = None
        for row in csv_reader:
            for num, s in enumerate(study_list):
                if s in row:
                    study_row_list[num].append(row[:-1])
    # make dictionaries for all studies
    for srows in study_row_list:
        rownum = len(srows)
        val_list_list = [[None for j in range(rownum)] for i in range(len(header_list))]
        for i, row in enumerate(srows):
            for j, val in enumerate(row):
                val_list_list[j][i] = val
        # make dictionary with header
        d = dict(zip(header_list, val_list_list))
        study_dict_list.append(d)

    return study_dict_list

def read_csv_headers(infile):
    """Return a list of csv headers as strings"""
    n = -1
    with open(infile,"r") as openfile:
        lines = openfile.readlines()
        for num, line in enumerate(lines):
            if "full_log" in line:
                n = num
            if num == n+1:
                hlist = line.split(',')[:-1]
    return hlist


def interpolate_gaps(vals, limit=None):
    """
    Fill gaps using linear interpolation, optionally only fill gaps up to
    a size of 'limit'

    """
    vals = np.asarray(vals)
    i = np.arange(vals.size)
    valid = np.isfinite(vals)
    filled = npinterp(i, i[valid], vals[valid])

    if limit is not None:
        invalid = ~valid
        for n in range(1, limit+1):
            invalid[:-n] &= invalid[n:]
        filled[invalid] = np.nan

    return filled

def fill_gaps(vals, dtype=None):
    """ Interpolate NAN in numpy array"""

    print("IN FIL GAP")
    if dtype == "time":
        for n, val in enumerate(vals):
            if val == 0:
                val = int((vals[n-1] + vals[n+1]) / 2)

    if dtype == None:
        for n, val in enumerate(vals):
            if val == 0:
                val = int((vals[n-1] + vals[n+1]) / 2)

    return 

if __name__ == "__main__":
    main()
