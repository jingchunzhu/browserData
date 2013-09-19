#!/usr/bin/env bash

for maf in $(ls ../../../PanCan/MAF/somatic_mafs_cleaned_filtered_syn1729383/*maf)
do
    filename=${maf##*/}
    tumor=${filename%.*}
    echo $maf, $tumor
    mafToSparseMatrix.py $maf | transpose.py > ../mutationsFromPanCancer/$tumor.txt
done
