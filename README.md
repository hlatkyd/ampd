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

Compile & Install
---
```
git clone http://www.github.com/hlatkyd/ampd.git
cd ampd
make all
```
Manually move the binaries to desired path.

Basic usage
---

```
ampd -f [infile] -o [outfile] -r [sampling_rate] -l [batch_length]
```
Returns the number of peaks to stdout.

### Options in detail
```
-h --help           Print help.
-v --verbose        Verbose output.
-f --infile         Input file, should only contain one float value each line
-o --outdir         Output directory. Depending on required output tpe, extensions
                    are added and multiple files are created. Also creates the 
                    required directories.
                    Default is [cwd]/ampd.out
-a --auxdir         Aux output directory containing all intermediate data. Should
                    be used for thorough manual checking only. Created in case
                    the option --output-all is set.
                    Default is [cwd]/ampd.aux
-r --sampling-rate  Sampling rate of the data in input file, in Hz.
                    Default is 100 Hz.
-l --batch_length   Large data files are processed in non-overlapping window approach
                    (batches). This parameter stands for the window length in seconds.
                    Default is 60s.
--preproc           Do preprocessing by applying simple filters. The utility program
                    ampdpreproc is called on the data with the specified filters.
--lpfilt            Lowpass filter in Hz.
--hpfilt            Highpass filter in Hz.
--output-rate       Output number of peaks per minute, for each batch window.
--output-peaks      Output peak indices corresponding to original data.
--output-meta       Output metadata to file.
--output-all        Output aux files, except local maxima scalogram.
--output-lms        Ouptut local maxima scalogram matrix in auxdir.
```
Option defaults are defined in ampd.h. Reset these as convenient, then recompile.

Utility programs
---
ampd needs an input file with a single value on each line, thus a little outside
peparation is needed. Some additional tools can be found in ./bin/ after compilation.

colextract, rowextract: prepare input file

ampdcheck: plot some outputs of ampd

### TODO
* clean up config file
* data type (eg.: respiration, ECG, pulsoxy) dependent defaults
* speeding up, MPI maybe?
* make ampd routine change windows adaptively in case peaks are incorrect
