#!/bin/bash

# 

datatype=$1
sr=100
testdatadir="/home/david/dev/ampd/test/data"
ampdroot="/home/david/dev/ampd"
testaux="$ampdroot/test/out/ampd.aux"
testout="$ampdroot/test/out/ampd.out"
echo ""
echo "Testing AMPD"
echo "---------------------------------------"
if [ -z "$1" ]; then
    echo "Argument needed. Try: 'resp', 'puls', 'robustness', 'resp_flip'"
    echo "See script code for more details"
fi
# Test respiration peak counting

if [ "$datatype" = "resp" ]; then 
    $ampdroot/bin/ampd -f $testdatadir/resp_raw.txt -v -l 60 \
        --output-all -r $sr --preproc -t resp -a $testaux -o $testout
    $ampdroot/scripts/ampdcheck.py $testaux/batch_0
fi
if [ "$datatype" = "resp_flip" ]; then 
    $ampdroot/bin/ampd -f $testdatadir/resp_flip_raw.txt -v -l 60 \
        --output-all -r $sr --preproc -t resp -a $testaux -o $testout
    $ampdroot/scripts/ampdcheck.py $testaux/batch_0
fi
# Test pulsoxy peak counting
if [ "$datatype" = "puls" ]; then 
    $ampdroot/bin/ampd -f $testdatadir/pulsoxy_raw.txt -v -l 60 \
        --output-all -r $sr --preproc -t puls -a $testoux -o $testout
    $ampdroot/scripts/ampdcheck.py $testaux/batch_0
fi
# test time and stability of peak count for different batch lengths
if [ "$datatype" = "robustness" ]; then 
    echo "testing pulsoxy data..."
    $ampdroot/bin/ampdpreproc -f $testdatadir/pulsoxy_raw.txt \
        -o $testdatadir/pulsoxy.txt -s $sr -h 2
    for i in 10 30 60 100
    do
        SECONDS=0
        printf "batch_length=$i, total peaks: "
        $ampdroot/bin/ampd -f $testdatadir/pulsoxy.txt -l $i -r $sr -o $testout
        printf "elapsed time: $SECONDS sec\n"
    done

    echo "testing resp data..."
    $ampdroot/bin/ampdpreproc -f $testdatadir/resp_raw.txt \
        -o $testdatadir/resp.txt -s $sr -h 0.2 -l 3
    for i in 10 30 60 100
    do
        SECONDS=0
        printf "batch_length=$i, total peaks: "
        $ampdroot/bin/ampd -f $testdatadir/resp.txt -l $i -r $sr -o $testout
        printf "elapsed time: $SECONDS sec\n"
    done
    echo "testing flipped resp data..."
    $ampdroot/bin/ampdpreproc -f $testdatadir/resp_flip_raw.txt \
        -o $testdatadir/resp_flip.txt -s $sr -h 0.2 -l 3
    for i in 10 30 60 100
    do
        SECONDS=0
        printf "batch_length=$i, total peaks: "
        $ampdroot/bin/ampd -f $testdatadir/resp_flip.txt -l $i -r $sr -o $testout
        printf "elapsed time: $SECONDS sec\n"
    done
fi
