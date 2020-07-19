#!/bin/bash

# 
echo ""
echo "Testing AMPD"
echo "---------------------------------------"
./bin/ampd -f test_data/resp.txt -v -l 120 --overlap=0.2 --output-all
#./bin/ampd -f test_data/pulsoxy.txt -v -l 10 --overlap=0.2 --output-all
#./util/ampdcheck ampd_out/batch_0
