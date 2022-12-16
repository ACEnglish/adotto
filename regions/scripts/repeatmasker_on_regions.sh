#!/bin/bash

~/scratch/misc_software/RepeatMasker/RepeatMasker -pa 16 -qq -e hmmer -species human -lcambig -nocut -div 50 -no_id \
    -s tr_regions.v0.3.fasta
