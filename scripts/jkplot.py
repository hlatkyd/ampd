#!/usr/bin/python3

import matplotlib.pyplot as plt
import numpy as np
import getopt
from jkfetch import fetch_study_log
import os
import sys

def usage():

    txt="""
    jkplot.py --start=[startid] --stop=[stopid] --exclude_file

    
    """

    print(txt)
    return 

# USER SETUP
#------------------------------------------------------------------------------
EXCLUDE_FILE = "exclude_id"
JK_PATH = "~/work/jk2020aug.csv"        # sequence log file
JK_HEADER_PATH = "~/work/csv_headers"   # header for sequence log

#------------------------------------------------------------------------------

optstr = "hv"
longopt = ["start=","stop=", "help", "exclude-file"]

def main():

    start = None
    stop = None
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
        elif o in ("--exclude-file"):
            exclude_list = read_exclude_file(EXCLUDE_FILE)
            #TODO
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
    if ((start != None) & (len(str(start)) not in [2,4,6,10,12])):
        error += 1 
    if ((stop != None) & (len(str(stop)) not in [2,4,6,10,12])):
        error += 1 
    if error > 1:
        print("ERROR: Wrong format for arguments --stop=arg or --start=arg, exiting..")
        print("start="+str(start)+" stop="+str(stop))
    if error > 0:
        usage()
        sys.exit(0)
    return

if __name__ == "__name__":
    main()
