#!/bin/bash

# 
echo ""
echo "Testing AMPD"
echo "---------------------------------------"
./bin/ampd -f test_data/extract_test.txt -v -a
#./util/ampdcheck ampd_out/batch_0
