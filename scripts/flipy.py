#!/usr/bin/python3
"""
Flip data along y axis.
Read from file, which contains one data values per line.
"""
import numpy as np
import os
import glob
import matplotlib.pyplot as plt
import argparse

def main():

    parser = argparse.ArgumentParser()
    parser.add_argument('path', type=str)
    args = parser.parse_args()

    if not os.path.isfile(args.path):
        print("Input is not a file!\n")
        exit()

    data = np.loadtxt(args.path, dtype=float)
    mean = np.mean(data)
    maxdata = np.amax(data)
    data = - data + maxdata
    out = "fliptest.txt"
    np.savetxt(out, data, fmt='%.3f')
    hist = np.histogram(data, bins=20)
    plt.hist(hist)
    plt.show()


        



if __name__ == "__main__":
    main()
