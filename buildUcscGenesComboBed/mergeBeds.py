#!/usr/bin/env python2.7

import argparse
from pycbio.hgdata import Bed

parser = argparse.ArgumentParser()
parser.add_argument('oldBed', type=str, help="Old input bed file")
parser.add_argument('newBed', type=str, help="New input bed file")
args = parser.parse_args()

#
# Make a dictionary of all transcripts (by name) from the old bed file
bedData = dict()
oldBedFp = open(args.oldBed)
for line in oldBedFp:
    line = line.rstrip()
    bb = Bed.Bed(line.split())
    bedData[bb.name] = bb
oldBedFp.close()

#
# Read the transcripts from the new bed file and store them in the dictionary.
# If any transcript exists in the old and new bed files, the data from the new
# bed file will replace the data from the old bed file.  This is what we want.
newBedFp = open(args.newBed)
for line in newBedFp:
    line = line.rstrip()
    bb = Bed.Bed(line.split())
    bedData[bb.name] = bb
newBedFp.close()

#
# Print out the contents of the dictionary, sorted by name.
for transcript in sorted(bedData.keys()):
    print bedData[transcript]
    

