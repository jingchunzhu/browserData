#!/usr/bin/env bash
dir=$(dirname ${BASH_SOURCE[0]})
/inside/home/jzhu/scripts/jre1.7.0_55/bin/java -jar -Xmx2480m $dir/../snpEff/snpEff.jar hg19 -c $dir/../snpEff/snpEff.config -quiet  -sequenceOntology -no-downstream -no-upstream -no-intergenic -no-intron -v -noStats $1 