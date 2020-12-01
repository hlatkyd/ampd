"""
plotutils.py

Utility functions used by pulsplot, respplot, logplot, etc
Keep in same directory with those scripts
"""

ID_PREFIX = "s_"    # prefix for study directory
ID_SUFFIX = ""      # suffix for study directory

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

#TODO delete
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
