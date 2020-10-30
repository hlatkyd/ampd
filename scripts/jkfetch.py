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
    Keys are found in jk header file 'csv_headers' under 'full_log'

    currently the keys are:

    ['studyid', 'pslabel', 'comment', 'scantime', 'tr', 'te', 'images', 'nt',
    'ss', 'seqid', 'fwhm', 'type', 'ratid', 'weight', 'measurement',
    'anesthesia', 'treatment', 'time', 'resp', 'bpm', 'isoflurane',
    'resp_var', 'bpm_var', 'op_expcomment', 'a_expcomment',
    'Prescan_FatOffset', 'H1offset', 'pwr90']

    """

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

    study_dict_list = []
    study_row_list = [[] for i in range(len(study_list))]
    header_list = read_csv_headers(header_path)
    print(header_list)

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

def getpath(a):

    if a[0] == "~":
        return os.path.expanduser('~')+a[1:]
    else:
        return os.path.abspath(a)


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
        if o in ("-i","--infile"):
            infile = getpath(a)
        elif o in ("-h","--header"):
            header = getpath(a)
        elif o == "-v":
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
