#!/usr/bin/env python

"""mafToSparseMatrix.py: given a MAF, generate a sparse matrix indicating
which genes have mutations of the selected type(s), as specified by
command line options (missense or nonsense by default"""

import argparse
import csv
import MySQLdb
import MySQLdb.cursors
from collections import OrderedDict
import re
import sys

def getHugoGeneSymbols():
    """Get a full list of current HUGO symbols"""
    symbols = list()
    db = MySQLdb.connect(read_default_file="~/.my.hg19.cnf")
    cursor = db.cursor(MySQLdb.cursors.DictCursor)
    cursor.execute("select symbol from hgnc order by symbol")
    for row in cursor.fetchall():
        symbol = row["symbol"]
        symbol = re.sub(".withdrawn$", "", symbol)
        if not (symbol in symbols):
            symbols.append(symbol)
    return(symbols)

def determineIfProteinCoding(listOfSymbols):
    """Given a list of gene symbols, return a dictionary indicating
    which ones code for protein."""
    isCodingGene = dict()
    db = MySQLdb.connect(read_default_file="~/.my.hg19.cnf")
    cursor = db.cursor(MySQLdb.cursors.DictCursor)
    for symbol in listOfSymbols:
        query = """select count(*) as codingTranscripts 
                    from knownGene kg, kgXref kx 
                    where kg.name = kx.kgID and kg.cdsStart != kg.cdsEnd 
                       and kx.geneSymbol = '%s'""" % (symbol)
        cursor.execute(query)
        row = cursor.fetchone()
        if row["codingTranscripts"] != 0:
            isCodingGene[symbol] = True
        else:
            isCodingGene[symbol] = False
    return(isCodingGene)


def printHeaderRow(symbols):
    """Print a header row with the HUGO symbols"""
    line = "Sample"
    for thisSymbol in symbols:
        line = line + "\t%s" % (thisSymbol)
    print line

def dataForNewSample(symbols):
    """Initialize a results vector for a new sample"""
    dataThisSample = OrderedDict()
    for thisSymbol in symbols:
        dataThisSample[thisSymbol] = 0
    return(dataThisSample)

def printDataThisSample(sampleId, symbols, dataThisSample):
    """Print a row of results for one sample"""
    line = sampleId
    for thisSymbol in symbols:
        line = line + "\t%d" % (dataThisSample[thisSymbol])
    print line

def main():
    mutationTypesIfCoding = "Nonsense|Missense|Nonstop|Frame_Shift|Splice_Site"
    mutationTypesIfNoncoding = mutationTypesIfCoding + "|In_Frame|RNA"
    parser = argparse.ArgumentParser()
    parser.add_argument("maf", type=str, help="Input MAF file")
    parser.add_argument("-d", "--debug", default=False, help="Debug?")
    args = parser.parse_args()

    hugoSymbols = getHugoGeneSymbols()
    isProteinCoding = determineIfProteinCoding(hugoSymbols)
    printHeaderRow(hugoSymbols)
    mafData = csv.DictReader(open(args.maf), delimiter="\t")
    currentSampleBarcode = ""
    for row in mafData:
        nextSampleBarcode = row["Tumor_Sample_Barcode"]
        if nextSampleBarcode != currentSampleBarcode:
            if currentSampleBarcode != "":
                printDataThisSample(currentSampleBarcode, hugoSymbols,
                                    sampleData)
            sampleData = dataForNewSample(hugoSymbols)
            currentSampleBarcode = nextSampleBarcode
        geneSymbol = row["Hugo_Symbol"]
        if (geneSymbol in hugoSymbols) and row["Mutation_Status"] == "Somatic":
            mutationType = row["Variant_Classification"]
            if isProteinCoding[geneSymbol]:
                if re.match(mutationTypesIfCoding, mutationType):
                    sampleData[geneSymbol] = 1
            else:
                if re.match(mutationTypesIfNoncoding, mutationType):
                    sampleData[geneSymbol] = 1
    printDataThisSample(currentSampleBarcode, hugoSymbols, sampleData)

if __name__ == "__main__":
    main()
        
