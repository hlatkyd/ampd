#!/usr/bin/python3

import os
import glob
import sys
import getopt

#------------------------------------------------------------------------------
#                               USER OPTIONS
#------------------------------------------------------------------------------

# data types to process present in input files: 
RESP = True
PULS = True
ECG = False

# should ampd generate aux output?
RESP_AUX = False
PULS_AUX = False
ECG_AUX = False
# sai data file format: [prefix][study_id][suffix]
PREFIX = "data_"
SUFFIX = ""
# study id format
#ID = s_$YEAR$MONTH$DAY$NUM, eg: s_2020040501

# colextract options
                # column numbering start at 0
RESP_COL = 1    # position of respiration data in raw sai data file, 
PULS_COL = 2    # position of pulsoxy data 
ECG_COL = 3     # ecg data column
# ampd options
LENGTH = 60     # batch-length
#saiproc defaults
VERBOSE_DEF = False

#------------------------------------------------------------------------------
optstr = "hv:i:o:"
optlong = ["help","verbose","edit","indir=","outdir="]

def usage():

    text="""
    Use ampd, rowextract, colextract to process MRI physiological data acquired by
    SA instruments system.

    Usage:
        saiproc.py -i [indir] -o [outdir]

    Optional arguments:
        -h --help
        -v --verbose

    Options:
        Various options should be set to accomodate input data format. See script
        source code and set these accordingly before calling.

    Edit in vim:
        saiproc.py --edit

    The imput directory should contain multi-column data files in text format. The
    format location of respiration, pulsoxy, etc timeseries and name of the files 
    should be specified within this script, as well as the parameters for ampd.
    The output layout:
    [outdir]/
            study_id_x/
                       resp.txt
                       puls.txt
                       resp.ampd.out/
                       puls.ampd.out/
                       resp.ampd.aux/       -- optional
                       puls.ampd.aux/       -- optional
                                batch_x
    """                         
    print(text)

def usage_simple():
    text="""
    Usage:
        saiproc.py -i [indir] -o [outdir]
    """

def params():
    """ Read parameters and create a dictionary"""

    par = dict()
    par["resp_col"] = RESP_COL
    par["puls_col"] = PULS_COL
    par["length"] = LENGTH
    par["verbose"] = VERBOSE_DEF
    par["prefix"] = PREFIX
    par["suffix"] = SUFFIX
    return par

def extract_study_id(full_path, par):
    """ Extract study id from the path to sai data file"""
    sid = os.path.basename(full_path)
    if par["prefix"] in sid and par["prefix"] != "":
        sid = sid.split(par["prefix"])[1]
    if par["suffix"] in sid and par["suffix"] != "":
        sid = sid.split(par["suffix"])[0]
    return sid

def which(program):
    """Imitates the unix 'which' command. Returns None if executable not in PATH"""
    def is_exe(fpath):
        return os.path.isfile(fpath) and os.access(fpath, os.X_OK)

    fpath, fname = os.path.split(program)
    if fpath:
        if is_exe(program):
            return program
    else:
        for path in os.environ["PATH"].split(os.pathsep):
            exe_file = os.path.join(path, program)
            if is_exe(exe_file):
                return exe_file

    return None

def main():

    par = params()
    if len(sys.argv) == 1:
        usage()
        sys.exit(0)
    try:
        opt, args = getopt.getopt(sys.argv[1:],optstr,optlong)
    except getopt.GetoptError as err:
        print(str(err))
        usage()
        sys.exit(2)
    outdir = None
    indir = None
    verbose = par["verbose"]
    for o, a in opt:
        if o in ("-v", "--verbose"):
            verbose = True
        elif o in ("-h", "--help"):
            usage()
            sys.exit(0)
        elif o == "--edit":
            if which("vim") is not None: 
                cmd = "vim "+str(os.path.realpath(__file__))
                os.system(cmd)
                sys.exit(0)
            else:
                print("Cannot find vim for editing. Open and edit the script manually.")
                sys.exit(0)
        elif o in ("-i", "--indir"):
            indir = str(a)
            indir = os.path.abspath(indir)
        elif o in ("-o", "--outdir"):
            outdir = str(a)
            outdir = os.path.abspath(outdir)
        else:
            assert False, "unhandled option"

    if indir == None or outdir == None:
        print("Both output and input directory arguments are required")
        usage_simple()
        sys.exit(2)
    # check for helper c routines
    if which("colextract") is None:
        print("ampd helper 'colextract' not found in PATH, exiting")
        sys.exit(2)
    if which("rowextract") is None:
        print("ampd helper 'rowextract' not found in PATH, exiting")
        sys.exit(2)

    # search for data files in input dir
    pattern = par["prefix"] + "*" + par["suffix"]
    file_list = sorted(glob.glob(indir+"/"+pattern))

    for num, f in enumerate(file_list):
        #TODO make batch processing easier to do??


        #TODO del, num>0 is only for testing
        print("Processing "+f+" "+str(num)+"/"+str(len(file_list)))
        """
        if num > 0:
            break
        """
        study_id = extract_study_id(f, par)
        outdir_study = outdir +"/"+ str(study_id)
        if os.path.isdir(outdir_study) == False:
            try:
                os.makedirs(outdir_study, exist_ok=True)
            except OSError as exc:
                if exc.errno == errno.EEXIST and os.path.isdir(outdir_study):
                    pass
                else:
                    raise

        # prepare ampd input
        resp_data = outdir_study + "/resp.txt"
        puls_data = outdir_study + "/puls.txt"
        ecg_data = outdir_study + "/ecg.txt"
        # resp
        cmd = "colextract"+" -f "+f +" -o "+resp_data + " -n "+str(par["resp_col"]) 
        os.system(cmd)
        # puls
        cmd = "colextract"+" -f "+f +" -o "+puls_data + " -n "+str(par["puls_col"]) 
        os.system(cmd)

        # prepare ampd output
        resp_ampd_out = outdir_study + "/resp.ampd.out"
        resp_ampd_aux = outdir_study + "/resp.ampd.aux"
        puls_ampd_out = outdir_study + "/puls.ampd.out"
        puls_ampd_aux = outdir_study + "/puls.ampd.aux"
        ecg_ampd_out = outdir_study + "/ecg.ampd.out"
        ecg_ampd_aux = outdir_study + "/ecg.ampd.aux"

        # run ampd
        if RESP:
            infile = resp_data
            cmd = "ampd"+" -f "+infile+" -o "+resp_ampd_out+" -t resp "+" -l 60"
            if RESP_AUX:
                cmd = cmd+" -a "+resp_ampd_aux +" --output-all"
            os.system(cmd)
        if PULS:
            infile = puls_data
            cmd = "ampd"+" -f "+infile+" -o "+puls_ampd_out+" -t puls "+" -l 60"
            if PULS_AUX:
                cmd = cmd+" -a "+puls_ampd_aux+" --output-all"
            os.system(cmd)
        if ECG:
            infile = ecg_data
            cmd = "ampd"+" -f "+infile+" -o "+ecg_ampd_out+" -t ecg "+" -l 60"
            if ECG_AUX:
                cmd = cmd+" -a "+ecg_ampd_aux+" --output-all"
            os.system(cmd)

    

    return


if __name__ == "__main__":
    main()
