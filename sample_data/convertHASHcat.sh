#!/bin/bash

# script extracts some of the columns from the HASH focal mechanism catalog

eqkfile=socal_1984_2003.hash.gz
ofile=${eqkfile%*.hash.gz}.llzsdrmi.txt

# get the info
zcat $eqkfile | \
    gawk '{lon=$9; lat=$8; z=$10; m=$11; str=$12; dip=$13; rake=$14; id=$7;
printf("%-10.4f %10.4f %6.3f %6.2f %6.1f %6.1f %6.1f %12i\n",
lon,lat,z,str,dip,rake,m,id);}' > $ofile

echo Written to $ofile
