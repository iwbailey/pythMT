#!/bin/bash

# Simple test to illustrate how psmecaSmtoSx.py may be used

psmecaSmtoSx=../psmecaSmtoSx.py

echo "Test help file... "
$psmecaSmtoSx -h
echo ""

function runtest()
{
    echo "In: $test_in"
    echo "Out: "$(echo $test_in | $psmecaSmtoSx)
    echo ""
}

echo "Test Conversion"
test_in="0 0 0 1 -1 0 0 0 0 20 0 0"
runtest

test_in="0 0 0 0 0 0 0 0 1 20 0 0"
runtest
