#!/bin/bash

# Simple test to illustrate how mtConvert may be used

mtconvert=../mtConvert.py

echo "Test help file... "
$mtconvert -h
echo ""

echo "Test SDR -> rotation axis"
s=0; d=45; r=-90
echo "Strike = $s, dip = $d, rake = $r"
echo "Output..."
echo 0 0 0 $s $d $r 1 | $mtconvert -i 2 -o 4 
