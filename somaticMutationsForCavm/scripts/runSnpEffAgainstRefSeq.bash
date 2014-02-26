#!/usr/bin/env bash
dir=$(dirname ${BASH_SOURCE[0]})
java -jar -Xmx2480m $dir/../snpEff/snpEff.jar hg19 -c $dir/../snpEff/snpEff.config -quiet  -sequenceOntolgy -v -noStats $1 