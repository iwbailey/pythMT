#!/bin/bash

# Check that psmeca format is read in correctly

head -5 ../sample_data/socal_cmts.psmeca | \
    python ../src/getMTinfo.py