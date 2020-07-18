#!/bin/bash

# 
echo ""
echo "Testing AMPD"
echo "---------------------------------------"
./bin/ampd -f test_data/resp.txt -v 2 -l 120 --overlap=0.2 --output-all
#./bin/ampd -f test_data/pulsoxy.txt -v -l 30 --overlap=0.2 --output-all
#./util/ampdcheck ampd_out/batch_0
