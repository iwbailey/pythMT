#!/bin/bash

# Test extracting some quantities using mtInfo.sh

mtinfo=../mtInfo.py

echo "Test help file... "
$mtinfo -h

echo -e "\nTest what happens when no arguments..."
$mtinfo

echo -e "\nTest reading from stdin..."
head -5 ../../sample_data/socal_cmts.psmeca | \
    $mtinfo -v

echo -e "\nTest reading from file..."
$mtinfo ../../sample_data/socal_cmts.psmeca -v

echo -e "\nTest --pos..."
head -5 ../../sample_data/socal_cmts.psmeca | \
    $mtinfo --pos

echo -e "\nTest --smt..."
head -5 ../../sample_data/socal_cmts.psmeca | \
    $mtinfo --smt

echo -e "\nTest --m0..."
head -5 ../../sample_data/socal_cmts.psmeca | \
    $mtinfo --m0

echo -e "\nTest --mw..."
head -5 ../../sample_data/socal_cmts.psmeca | \
    $mtinfo --mw

echo -e "\nTest --lognm..."
head -5 ../../sample_data/socal_cmts.psmeca | \
    $mtinfo --lognm

echo -e "\nTest --fclvd..."
head -5 ../../sample_data/socal_cmts.psmeca | \
    $mtinfo --fclvd

echo -e "\nTest --Gamma..."
head -5 ../../sample_data/socal_cmts.psmeca | \
    $mtinfo --Gamma

echo -e "\nTest --frr..."
head -5 ../../sample_data/socal_cmts.psmeca | \
    $mtinfo --frr

echo -e "\nTest several things together"
head -5 ../../sample_data/socal_cmts.psmeca | \
    $mtinfo --pos --smt --mw --fclvd

echo -e "\nTest paste"
head -5 ../../sample_data/socal_cmts.psmeca | \
    $mtinfo --paste --fclvd

exit
