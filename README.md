Automatic Multi-scale Peak Detection
====
A robust algorith for peak detection implemented in c. Intended usage is
finding and counting peaks and peak rates in rat physiological datasets, such
as respiration and pulsoxymetry measured by SA Instruments Small Animal
Monitoring System during an MRI session.

Reference paper:
Scholkmann et. al:
An Efficient Algorithm for Automatic Peak Detection in Noisy Periodic and
Quasi-Periodic Signals
doi:10.3390/a5040588

<img src=https://raw.githubusercontent.com/hlatkyd/ampd/master/doc/fig1.png"
alt="example respiration peak detection" width="1200"/>

Compile & Install
---
```
git clone http://www.github.com/hlatkyd/ampd.git
cd ampd
make
make install
```
Manually move the helper binaries and scripts to desired path if needed.

Basic usage
---

```
ampd -f [infile] -o [outdir] -r [sampling_rate] -l [batch_length]
```
Returns the number of peaks to stdout. It is recommended to specify ```--datatype``` 
in production and leave ```--sampling-rate``` and ```batch-length``` as defaults in
most cases, when the input data in question is respiration or pulsoxymetry coming 
from the SA Instuments system.
```
ampd -f [infile] -o [outdir] -t [datatype]
```

### Options in detail
```
-h --help           Print help.
-v --verbose        Verbose output.
-f --infile         Input file, should only contain one float value each line
-o --outdir         Output directory. Depending on required output tpe, extensions
                    are added and multiple files are created. Also creates the 
                    required directories.
                    Default is [cwd]/ampd.out
-t --datatype       Option to preset preprocessing defaults specified for a datatpye,
                    Right now either "resp" or "puls" is accepted, otherwise it is
                    ignored. Calling with --lpfilt=x, --hpfilt=x, etc overrides these
                    defaults.
-a --auxdir         Aux output directory containing all intermediate data. Should
                    be used for thorough manual checking only. Created in case
                    the option --output-all is set.
                    Default is [cwd]/ampd.aux
-r --sampling-rate  Sampling rate of the data in input file, in Hz.
                    Default is 100 Hz.
-l --batch-length   Large data files are processed in non-overlapping window approach
                    (batches). This parameter stands for the window length in seconds.
                    Default is 60s.
--preproc           Do preprocessing by applying simple filters.
--lpfilt            Lowpass filter in Hz.
--hpfilt            Highpass filter in Hz.
--output-rate       Output number of peaks per minute, for each batch window.
--output-peaks      Output peak indices corresponding to original data.
--output-meta       Output metadata to file.
--output-all        Output aux files, except local maxima scalogram.
--output-lms        Ouptut local maxima scalogram matrix in auxdir.
--autoflip          Flip data along Y axis if events are minima as determined by
                    the centre of mass of a histogram of the data. Default is ON
```
Some defaults in case optional arguments are not given are defined in ampd.h.
Reset these as convenient, then recompile.


###Testing

Test data samples are found in ```./test/data``` directory. The script ```test.sh``` can be called with arguments corrspondind to the type of data, currently: "resp" and "puls". See code for more detail.

Utility programs
---
ampd needs an input file with a single value on each line, thus a little outside
peparation is needed. Some additional tools can be found in ./bin/ after compilation.

colextract, rowextract: prepare input file

ampdcheck.py:   plot some outputs of ampd
flipy.py:        flip data along Y axis

### TODO
* clean up config file
* data type (eg.: respiration, ECG, pulsoxy) dependent defaults
* speeding up, MPI maybe?
* make ampd routine change windows adaptively in case peaks are incorrect
