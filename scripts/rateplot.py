#!/usr/bin/python3

import matplotlib.pyplot as plt
from matplotlib.ticker import FormatStrFormatter
import numpy as np
import getopt
import sys
import os
import glob
from pathlib import Path    # Python > 3.5 
import csv

# User setup; modify as needed
#--------------------------------------------------------------
# datatype to plot. Files are searched wirh names DTYPE[i].rate
DTYPE = ["resp","puls"]
DTYPE_THRESH = [(30,150),(100,500)]
AUTO_EXCLUDE = 1     # exlude where data is not within bounds
AUTO_EXCLUDE_THRESH = 0.01 # acceptable ratio of data out of bounds
# study directories are searched as $IS_PREFIX[YYYYMMDDNN]$ID_SUFFIX
ID_PREFIX = "s_"    # prefix for study directory
ID_SUFFIX = ""      # suffix for study directory
NUM = len(DTYPE)    # number of data type, plots
HEADER_LENGTH = 2   # maximum number of header lines in .rate files

#--------------------------------------------------------------
optstr = "hi:v"
longopt = ["start=","stop=","indir=", "help", "exclude="]


def usage():
    text="""
    Usage:
        rateplot.py -i [indir] --start=[start_date] --stop=[stop_date] 

        Optional arguments:
        --exclude=[date1, date2, ...]
        --pattern=[rate_file_pattern]   // unused
        -v --verbose
        -h --help

    Plot multiple respiration, pulse, etc rates on each other.

    Input directory should contain the study directories named 's_[YYYYMMDDNN]',
    eg: s_2020050201. Study directories are searched for ampd output files with
    the suffix 'rate', eg: resp.rate. These are the main source of the plots.

    The optional arguments start and stop should have the following syntax:
    - 2 numbers: (eg '--start=01 --stop=06') means the ID on the last day
    - 4 numbers: (eg: '--start=0201' --stop=0405) means the studies in the last
      month from day 02 ID 01 to day 04 ID 05
    - 6 numbers: represent months as well in the last year
    - 10 numbers: represent the full IDs
    
    The optional argument '--exclude' uses these same formats. Multiple files
    can be excluded using comma as separator, eg: --exclude=073001,062002, ...
    """
    print(text)

def check_date_format(date):
    """Return True if start, stop, exclude arguments are of the correct format"""
    if date == None:
        return False
    if len(str(date)) not in [2,4,6,10]:
        return False
    if not isdigit(date):
        return False
    return True

def clean_arg_list(arg_orig):
    """ Return list of exclude dates"""
    arg = arg_orig
    # clean string
    chars = [" ","\t","{","}","[","]","(",")"]
    for c in chars:
        arg.replace(c,"")
    # arrange in list
    if ',' in arg:
        ret = arg.split(',')
    else:
        ret = [arg]
    # check if digits
    for item in ret:
        if not item.isdigit():
            print("ERROR: wrong arument syntax: "+str(arg_orig))
            sys.exit()
    return ret

def line_hover(event):
    """ interactive highlight"""
    ax = plt.gca()
    for line in ax.get_lines():
        if line.contains(event):
            pass

def generate_exclude_text(slist):

    text = "Exluded\n"
    for s in sorted(slist):
        text += os.path.basename(s)+"\n"
    return text

def main():

    start = None
    stop = None
    indir = None
    output = None
    verbose = False
    exclude_list = []
    try:
        opts, args = getopt.getopt(sys.argv[1:], optstr, longopt)
    except getopt.GetoptError as err:
        # print help information and exit:
        print(str(err)) 
        usage()
        sys.exit(2)
    for o, a in opts:
        if o == "-v":
            verbose = True
        elif o in ("-h", "--help"):
            usage()
            sys.exit()
        elif o in ("-i", "--indir"):
            indir = a
        elif o in ("--start"):
            start = a
        elif o in ("--stop"):
            stop = a
        elif o in ("--exclude"):
            exclude_list = clean_arg_list(a)
            print("exclude="+str(exclude_list))
        else:
            assert False, "unhandled option"

    # check for correct input arguments
    error = 0
    if len(sys.argv) == 1:
        usage()
        sys.exit(0)

    if indir == None:
        print("ERROR: Please specify input directory.")
        error += 1
    elif not os.path.isdir(indir):
        print("ERROR: Specified input directory does not exist or is not a "+
                "directory.Exiting...")
        sys.exit(0)
    if ((start != None) & (len(str(start)) not in [2,4,6,10])):
        error += 1 
    if ((stop != None) & (len(str(stop)) not in [2,4,6,10])):
        error += 1 
    if error > 1:
        print("ERROR: Wrong format for arguments --stop=arg or --start=arg, exiting..")
        print("start="+str(start)+" stop="+str(stop))
    if error > 0:
        usage()
        sys.exit(0)


    #---------- gather files to plot-------------------------------------------
    pattern = ID_PREFIX + "*"+ID_SUFFIX
    study_dir_list = []
    exclude_study_list = []
    data_files = [[] for i in range( NUM )]
    count = 1 if stop == None else 0

    for studydir in sorted(glob.glob(indir+"/"+pattern), reverse=True):

        if count == 0:
            if studydir.endswith(stop): # start counting back from stop
                count = 1
            else:
                continue
        # check for excluded
        if studydir.endswith(tuple(exclude_list)):
            exclude_study_list.append(studydir)
            continue
        study_dir_list.append(studydir)
        # search files
        for path in Path(studydir).rglob('*.rate'):
            for i in range(NUM):
                if DTYPE[i] in str(path.name):
                    data_files[i].append(path)

        if start is not None:
            if studydir.endswith(start):# stop counting once 'start' is encountered
                break

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

    # ---------------------------auto excldue --------------------------------
    # 1 if indexed data is to be plotted, change to 0 if excluded
    plot_list = [1 for i in range(len(study_dir_list))]
    for num, study, in enumerate(study_dir_list):
        for i in range(len(DTYPE)):

            d = data_arr[i][num]
            d = d[ d < DTYPE_THRESH[i][1]]
            d = d[ d > DTYPE_THRESH[i][0]]
            if not len(d) > len(data_arr[i][num]) * AUTO_EXCLUDE_THRESH:
                plot_list[num] = 0
                print("Exlcuding "+str(study_dir_list[num]))


    # --------------------------------plotting---------------------------------
    
    # check for samebatch lengts and sampling rates
    bl = batch_lengths[0]
    sr = sampling_rates[0]
    for i in range(len(batch_lengths)):
        if batch_lengths[i] != bl:
            print("WARNING: Batch lengths are not the same.")
        if sampling_rates[i] != sr:
            print("WARNING: Sampling rates are not the same.")
    study_list = [os.path.basename(x) for x in study_dir_list]
    
    fig, axes = plt.subplots(NUM,1, figsize=(15,7))
    data_xmax = 0
    for num, ax in enumerate(axes):
        for n, data in enumerate(data_arr[num]):
            ax.plot(data, label=study_list[n])
            if len(data) > data_xmax:
                data_xmax = len(data)
        maj_tick_freq = batch_length
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
            text = generate_exclude_text(exclude_study_list)
            props = dict(boxstyle='square', facecolor='white', alpha=0.5)
            ax.legend(loc='upper left',bbox_to_anchor=(1.,1))
            ax.text(1.2, 0.95, text, verticalalignment='top',fontsize=10,
                    bbox=props,transform=ax.transAxes)
    plt.show()
            
if __name__ == "__main__":
    main()
