#!/usr/bin/python3
"""
Flip data along y axis.
Read from file, which contains one data values per line.
"""
import numpy as np
import os
import glob
import getopt
import sys

def usage():
    txt = """
    flipy.py
    --------
    Flip data by applying transform to each points:
    data = - data + 2 * mean
    Input file is a text file containing float in strictly one column.

    Usage:
    flipy.py -f [infile] -o [outfile]
    """
    print(txt)

def main():

    if len(sys.argv) == 1:
        usage()
        sys.exit()
    optstr = "ho:f:"
    longopt = ["help", "infile=","outfile="]
    try:
        opts, args = getopt.getopt(sys.argv[1:], optstr, longopt)
    except getopt.GetoptError as err:
        print(str(err))
        usage()
        sys.exit()
    for o, a in opts:
        if o in ("-h", "--help"):
            usage()
            sys.exit()
        if o in ("-f", "--infile"):
            infile = str(a)
        if o in ("-o", "--outfile"):
            outfile = str(a)

    if not os.path.isfile(infile):
        print("Input is not a file!\n")
        usage()
        sys.exit()
    if not os.path.isabs(outfile):
        outfile = os.getcwd() + "/" + outfile

    data = np.loadtxt(infile, dtype=float)
    mean = np.mean(data)
    data = - data + 2 * mean
    np.savetxt(outfile, data, fmt='%.3f')
    #hist = np.histogram(data, bins=20)
    #plt.hist(hist)
    #plt.show()

if __name__ == "__main__":
    main()
