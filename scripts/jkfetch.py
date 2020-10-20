#!/usr/bin/python3

import numpy as np
import os
import csv
import getopt
import sys
"""
Functions to grab data from jk log files
"""

def fetch_study_log(study_list, jk_path, header_path):
    """
    Return a list a dictionaries, one dictionary per study.
    One dictionary contains lists of various parameters, one element 
    corresponding to a sequence.
    Keys:
    'seq'

    """

    study_dict_list = []
    for s in study_list:
        s_dict = {"name":s}
        study_dict_list.append(s_dict)

    with open(jk_path) as csv_file:
        csv_reader = csv.reader(csv_file,delimiter=',')
        line_count = 0
        for row in csv_reader:
            # 1 row represents 1 sequence in jk.csv
            if any(s in row for s in sorted(study_list,reverse=True)):
                print(s)



    return study_dict_list

def getpath(a):

    if a[0] == "~":
        return os.path.expanduser('~')+a[1:]
    else:
        return os.path.abspath(a)

def read_csv_headers(infile):
    """Return a list of csv headers as strings"""

    hlist = []
    with open(infile,) as openfile:


    return hlist
    

def main():

    study_list_test = ["s_2020070701", "s_2020070702"]
    infile = None
    header = None
    verbose = False
    try:
        opts, args = getopt.getopt(sys.argv[1:],"i:h:v",["infile=","header="])
    except getopt.GetoptError as err:
        print(str(err))
        sys.exit(2)
    for o, a in opts:
        if o in ("i","--infile"):
            infile = getpath(a)
        elif o in ("h","--header"):
            header = getpath(a)
        elif o == "v":
            verbose = True
        else:
            assert False, "unhandled option"
    if (infile == None or header == None):
        print("wrong input")
        sys.exit(2)

    fetch_study_log(study_list_test, infile, header)

    return

if __name__ == "__main__":

    main()
