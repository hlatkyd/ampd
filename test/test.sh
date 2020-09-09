#!/bin/bash

# 

datatype=$1
sr=100
echo ""
echo "Testing AMPD"
echo "---------------------------------------"
if [ -z "$1" ]; then
    echo "Argument needed. Try: 'resp', 'puls', 'robustness'"
fi


if [ "$datatype" = "resp" ]; then 
    #./bin/ampdpreproc -v -f test_data/resp_raw.txt -o test_data/resp.txt -s $sr -h 0.2 -l 3
    ./bin/ampd -f test_data/resp_raw.txt -v -l 60 --overlap=0.0 --output-all -r $sr --preproc -t resp
    ./scripts/ampdcheck.py ampd.aux/batch_0
fi
#./bin/ampd -f test_data/pulsoxy.txt -v -l 10 --overlap=0.2 --output-all
#./util/ampdcheck ampd_out/batch_0
if [ "$datatype" = "puls" ]; then 
    #./bin/ampdpreproc -v -f test_data/pulsoxy_raw.txt -o test_data/pulsoxy.txt -s $sr -h 2 
    ./bin/ampd -f test_data/pulsoxy_raw.txt -v -l 60 --overlap=0.0 --output-all -r $sr --preproc -t puls
    ./scripts/ampdcheck.py ampd.aux/batch_0
fi
# test time and stability of peak count for different batch lengths
if [ "$datatype" = "robustness" ]; then 
    echo "testing pulsoxy data..."
    ./bin/ampdpreproc -f test_data/pulsoxy_raw.txt -o test_data/pulsoxy.txt -s $sr -h 2
    for i in 10 30 60 100
    do
        SECONDS=0
        printf "batch_length=$i, total peaks: "
        ./bin/ampd -f test_data/pulsoxy.txt -l $i -r $sr
        printf "elapsed time: $SECONDS\n"
    done

    echo "testing resp data..."
    ./bin/ampdpreproc -f test_data/resp_raw.txt -o test_data/resp.txt -s $sr -h 0.2 -l 3
    for i in 10 30 60 100
    do
        SECONDS=0
        printf "batch_length=$i, total peaks: "
        ./bin/ampd -f test_data/resp.txt -l $i -r $sr
        printf "elapsed time: $SECONDS\n"
    done
fi