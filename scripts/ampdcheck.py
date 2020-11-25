#!/usr/bin/python3
"""
ampdcheck

Utility to plot the results of AMPD peak finding algorithm along with the
input and intermediate steps as well. Input can either be a singular file
or the directory cntaining all ampd aux out files.

Usage #1:
    ampdcheck [path/to/aux_dir/batch_dir]

    This aux output of ampd consist of the files:
        raw.dat
        detrend.dat
        lms.dat
        gamma.dat
        sigma.dat
        peaks.dat
        peaknum

    ampdcheck reads these files and plots all the data in orderly manner

Usage #2:
    ampdcheck [path/to/aux_dir]

Usage #3:
    ampdcheck [path/to/file]

    Simply plot the data as vector or matrix, whatever it finds. 

"""
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.widgets import Slider
import argparse
import csv
import os
import glob

DATFILES = ["raw.dat","detrend.dat","gamma.dat","sigma.dat",\
            "peaks.dat","param.txt","smoothed.dat"]
# linewidth of timeseries and gamma
LW_DEF = 1
# alpha of peak lines
A_DEF = 0.5

def main():

    parser = argparse.ArgumentParser()
    parser.add_argument('path', type=str, help='input dat file')
    #parser.add_argument('--lms', dest='lms',action='store_const',\
    #                    default=0,const=1,nargs=0,help='option to plot LMS')
    args = parser.parse_args()
    # check input path
    args.path = _abspath(args.path)

    # if argument is only a single file, plot it and return
    if os.path.isfile(args.path):
        check_single_input(args.path)
        return 0

    elif _is_batch_path(args.path):
        batches = sorted(glob.glob(os.path.dirname(args.path)+"/batch*"))
        n_batches = len(batches) 
        start_batch_num = int(os.path.basename(args.path).split("batch_")[-1] )
        start_batch_path = args.path

    elif _is_aux_path(args.path):
        batches = sorted(glob.glob(args.path+"/batch*"))
        n_batches = len(batches) 
        start_batch_path = batches[0]
        start_batch_num = 0


    # prepare figure
    fig, ax = batch_plot(start_batch_path)

    # add slider
    axis_color = "lightgoldenrodyellow"
    slider_init = 5
    slider_ax = fig.add_axes([0.2, 0.01, 0.60, 0.03], facecolor=axis_color)
    slider = Slider(slider_ax, 'batch',0, n_batches, valinit=start_batch_num,valstep=1)
    def slider_on_changed(mouse_event):
        fig.canvas.draw_idle()
        new_batch_path = _get_batch(start_batch_path, int(slider.val))
        batch_update(fig, ax, new_batch_path)
    slider.on_changed(slider_on_changed)
    # plot figure
    plt.show()

    return 0


def batch_update(fig, ax, batch_path):
    """Update axes from data in another batch"""
    pdict = load_param(batch_path+"/param.txt")
    for f in DATFILES:
        full_path = batch_path + "/" + f
        if f != "param.txt":
            data = np.loadtxt(full_path, delimiter='\t')
            x_axis_max = len(data) / int(float(pdict["sampling_rate"]))
            x_axis = np.arange(0,x_axis_max, len(data)*float(pdict["sampling_rate"]))
        if f == "raw.dat":
            ax[0].clear()
            ax[0].plot(data,linewidth=LW_DEF)
            datanum = len(data)
        elif f == "smoothed.dat":
            ax[1].clear()
            ax[1].plot(data,linewidth=LW_DEF)
        elif f == "detrend.dat":
            ax[2].clear()
            ax[2].plot(data,linewidth=LW_DEF)
        elif f == "gamma.dat":
            ax[3].clear()
            ax[3].plot(data,linewidth=LW_DEF)
        elif f == "sigma.dat":
            ax[4].clear()
            ax[4].plot(data,linewidth=LW_DEF)
        elif f == "peaks.dat":
            if data.ndim == 0: # hack if only one peak was found
                data = [data]
            for point in data:
                ax[0].axvline(x=point,color='r',zorder=0,alpha=A_DEF)
                #print("ylim = " + str(ax[0].get_ylim()))
                #ax[1].axvline(x=point,color='r',zorder=0)
                ax[2].axvline(x=point,color='r',zorder=0,alpha=A_DEF)
            n_peaks = len(data)
        else:
            pass

    # format axes
    for i in range(5):
        if i is not 4:
            ax[i].get_xaxis().set_ticks([])
        ax[i].margins(x=0)
    #plot rest based on params
    # plot vertical line on gamma min = lambda
    ax[3].axvline(int(pdict["lambda"]),color="orange")
    sampling_rate = float(pdict["sampling_rate"])
    x = np.linspace(0,datanum, datanum)
    y = float(pdict["fit_a"]) * x/sampling_rate + float(pdict["fit_b"])
    ax[0].plot(x, y, ':r')

    """
    # util formatting
    text = ['raw','detrend','smoothed','gamma','sigma']
    for i in range(n):
        ax[i].text(0.5,0.87,text[i],horizontalalignment='center',
                transform=ax[i].transAxes)
        #ax[i].set_xticklabels([])
    """
    # finishing up figure
    sampling_rate = 100
    lambdaa = pdict["lambda"]
    fig.suptitle(batch_path+"\n"+"sampling_rate=" + str(sampling_rate)\
                +", lambda=" + str(lambdaa) + ", n_peaks="+str(n_peaks))
    fig.canvas.draw()

def batch_plot(path):
    """
    Plot raw, smoothed, gamma, sigma, peaks first in one fig,
    the plot LMS in another fig
    """
    # the usual output of ampd, these should exist within the batch dir

    # check if all other exists
    for f in DATFILES:
        if f not in os.listdir(path):
            print("Cannot find file '"+f+"' exiting...\n")
            quit()
    # figure setup 
    n = len(DATFILES) - 2
    fig, ax= plt.subplots(n, 1, figsize=(14,7))

    # load param dictionary
    pdict = load_param(path+"/param.txt")

    plot_list = [None] * n
    for f in DATFILES:
        full_path = path + "/" + f
        if f != "param.txt":
            data = np.loadtxt(full_path, delimiter='\t')
            x_axis_max = len(data) / int(float(pdict["sampling_rate"]))
            x_axis = np.arange(0,x_axis_max, len(data)*float(pdict["sampling_rate"]))
        if f == "raw.dat":
            p, = ax[0].plot(data,linewidth=LW_DEF)
            datanum = len(data)
        elif f == "smoothed.dat":
            p, = ax[1].plot(data,linewidth=LW_DEF)
        elif f == "detrend.dat":
            p, = ax[2].plot(data,linewidth=LW_DEF)
        elif f == "gamma.dat":
            p, = ax[3].plot(data,linewidth=LW_DEF)
        elif f == "sigma.dat":
            p, = ax[4].plot(data,linewidth=LW_DEF)
        elif f == "peaks.dat":
            if data.ndim == 0: # hack if only one peak was found
                data = [data]
            for point in data:
                ax[0].axvline(x=point,color='r',zorder=0,alpha=A_DEF)
                #print("ylim = " + str(ax[0].get_ylim()))
                #ax[1].axvline(x=point,color='r',zorder=0)
                ax[2].axvline(x=point,color='r',zorder=0,alpha=A_DEF)
            n_peaks = len(data)
        else:
            pass

    # format axes
    for i in range(5):
        if i is not 4:
            ax[i].get_xaxis().set_ticks([])
        ax[i].margins(x=0)

    #plot rest based on params
    # plot vertical line on gamma min = lambda
    ax[3].axvline(int(pdict["lambda"]),color="orange")

    sampling_rate = float(pdict["sampling_rate"])
    x = np.linspace(0,datanum, datanum)
    y = float(pdict["fit_a"]) * x/sampling_rate + float(pdict["fit_b"])
    ax[0].plot(x, y, ':r')

    # util formatting
    text = ['raw','detrend','smoothed','gamma','sigma']
    for i in range(n):
        ax[i].text(0.5,0.87,text[i],horizontalalignment='center',
                transform=ax[i].transAxes)
        #ax[i].set_xticklabels([])

    plt.subplots_adjust(wspace=0,hspace=0)

    """
    Plot LMS
    if exists_lms == 1:
        fig2, ax2 = plt.subplots(1,2,figsize=(10,7))
        # plot lms
        full_path = path + "/" + "lms.dat"
        data = np.loadtxt(full_path, delimiter='\t')
        ax2[0].imshow(data,aspect='auto')
        _lambda = int(pdict["lambda"])
        ax2[1].imshow(data[:_lambda,:],aspect='auto')
        plt.tight_layout()


    """
    # finishing up figure
    sampling_rate = 100
    lambdaa = pdict["lambda"]
    fig.suptitle(path+"\n"+"sampling_rate=" + str(sampling_rate)\
                +", lambda=" + str(lambdaa) + ", n_peaks="+str(n_peaks))

    return fig, ax

def load_param(paramfile):
    """
    Return a dictionary from param.txt, given as input which contains
    ampd_param struct memebers in a format name=val at each line.

    """
    name = ["sampling_rate", "datatype", "a", "rnd_factor", \
            "fit_a", "fit_b", "fit_r", "lambda", "sigma_thresh",\
            "peak_thresh"]
    val = [None] * len(name)
    with open(paramfile, "r") as f:
        lines = f.read().splitlines()
        for line in lines:
            sline = line.split('=')
            for j in range(len(name)):
                if name[j] == sline[0]:
                    if len(sline) == 2:
                        val[j] = sline[1]
                    elif len(sline) == 1:
                        val[j] = None

    return dict(zip(name, val))


def _is_aux_path(path):
    """Return True if path is ampd aux directory, containing the batch direcotries"""
    contents = sorted(glob.glob(path+"/batch*"))
    if len(contents) != 0:
        return True

    return False

def _is_batch_path(path):
    """Return true if input path is a batch directory"""

    # should contain specific filenames
    file_list = ["raw.dat","smoothed.dat","detrend.dat","param.txt","sigma.dat",\
                "peaks.dat"]
    for f in file_list:
        if not os.path.isfile(path+"/"+f):
            return False

    return True

def check_single_input(path):
    """
    File argument only as input, plot it
    """
    if os.path.isfile(path):
        data = np.loadtxt(path, delimiter='\t')
        if len(data.shape) == 1:
            plt.plot(data)
            plt.margins(x=0)
            plt.tight_layout()
            plt.show()
        elif len(data.shape) == 2:
            plt.imshow(data)
            plt.tight_layout()
            plt.show()
        else:
            print("Cannot display data with shape "+str(data.shape))
    else:
        return 0

def _get_batch(start, val):
    """Return the path to batch directory specified by integer 'val'. """
    return os.path.dirname(start) + "/batch_"+str(val)

def _abspath(path):
    """Return absolute path. Handle ~ for home dir. """
    if path[0] == '~':
        return os.path.expanduser(path)
    else:
        return os.path.abspath(path)

if __name__ == "__main__":
    main()

