#!/bin/bash

# 
echo ""
echo "Testing AMPD"
echo "---------------------------------------"
#./bin/ampd -f test_data/resp.txt -v -a -l 120 -p 0.2
./bin/ampd -f test_data/pulsoxy.txt -va -l 30 -p 0.2
#./util/ampdcheck ampd_out/batch_0
