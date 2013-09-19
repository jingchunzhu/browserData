#!/usr/bin/env python2.7

import argparse
from pycbio.hgdata import Bed

parser = argparse.ArgumentParser()
parser.add_argument('--xref', nargs='+',
                    help="list of all UCSC Genes xref files, oldest first")
args = parser.parse_args()

#
# Make a dictionary to mapp transcript accessions to gene symbols 
accessionToSymbol = dict()

#
# Read transcript accession and gene symbol pairs from the xref file and store
# them in the dictionary.  Read the files in command line argument order, which
# should be oldest to newest.  In cases where the gene symbol for an accession
# has changed, this will have the effect of superceding the older mapping with
# the newer one.  Since the newer mappins are more accurate, this is what we
# want.
for thisXrefFile in args.xref:
    thisXrefFp = open(thisXrefFile)
    for line in thisXrefFp:
        tokens = line.rstrip().split('\t')
        accession = tokens[0]
        geneSymbol = tokens[4]
        accessionToSymbol[accession] = geneSymbol
    thisXrefFp.close()


#
# Print out the contents of the dictionary, sorted by transcript accession.
for accession in sorted(accessionToSymbol.keys()):
    print accession, accessionToSymbol[accession]
    

