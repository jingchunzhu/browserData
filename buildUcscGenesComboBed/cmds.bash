#!/usr/bin/env bash

#
# Build a hg18 and an hg19 bed file that contains all UCSC Genes accessions
# that have ever been.
mergeBeds.py /hive/data/genomes/hg19/bed/ucsc.11/ucscGenes.bed \
 /hive/data/genomes/hg19/bed/ucsc.12/ucscGenes.bed > ucscGenes.hg19.11_12.bed
mergeBeds.py ucscGenes.hg19.11_12.bed \
 /hive/data/genomes/hg19/bed/ucsc.13/ucscGenes.bed \
 > ucscGenes.hg19.11_12_13.bed
liftOver ucscGenes.hg19.11_12_13.bed \
  /hive/data/genomes/hg19/bed/liftOver/hg19ToHg18.over.chain.gz \
  ucscGenes.hg19.hg18_lifted.11_12_13.bed /dev/null

mergeBeds.py /hive/data/genomes/hg18/bed/ucsc.10/ucscGenes.bed \
  /hive/data/genomes/hg18/bed/ucsc.11/ucscGenes.bed > ucscGenes.hg18.10_11.bed
liftOver ucscGenes.hg18.10_11.bed \
  /hive/data/genomes/hg18/bed/liftOver/hg18ToHg19.over.chain.gz \
  ucscGenes.hg18.hg19_lifted.10_11.bed /dev/null

mergeBeds.py ucscGenes.hg18.hg19_lifted.10_11.bed ucscGenes.hg19.11_12_13.bed \
  > ucscGenes.hg19.bed
mergeBeds.py ucscGenes.hg19.hg18_lifted.11_12_13.bed ucscGenes.hg18.10_11.bed \
  > ucscGenes.hg18.bed

# 
# Build a mapping of accession name to gene symbol
(\
  cat /hive/data/genomes/hg19/bed/ucsc.13/ucscGenes.xref \
  | awk -F'\t' '{ print $1 "\t" $5}' | sort ; \
  cat /hive/data/genomes/hg19/bed/ucsc.12/ucscGenes.xref \
  | awk -F'\t' '{ print $1 "\t" $5}' | sort ; \
  cat /hive/data/genomes/hg19/bed/ucsc.11/ucscGenes.xref \
  | awk -F'\t' '{ print $1 "\t" $5}' | sort ; \
  cat /hive/data/genomes/hg18/bed/ucsc.11/ucscGenes.xref \
  | awk -F'\t' '{ print $1 "\t" $5}' | sort ; \
  cat /hive/data/genomes/hg18/bed/ucsc.10/ucscGenes.xref \
  | awk -F'\t' '{ print $1 "\t" $5}' | sort  \
) |sort |uniq > ucscGenes.accession.symbol.tab 

#
# For the probe map, extract the name, chrom, start, end, and strand from
# the bed files.  Sort them, and join with the accession to symbol file.
cat ucscGenes.hg19.bed |awk '{ print $4 "\t" $1 "\t" $2 "\t" $3 "\t" $6 }' |sort \
  | join -o 1.1,1.2,2.2,2.3,2.4,2.5 ucscGenes.accession.symbol.tab - \
> ucscGenesHg19Probe
cat ucscGenes.hg18.bed |awk '{ print $4 "\t" $1 "\t" $2 "\t" $3 "\t" $6 }' |sort \
  | join -o 1.1,1.2,2.2,2.3,2.4,2.5 ucscGenes.accession.symbol.tab - \
> ucscGenesHg18Probe
