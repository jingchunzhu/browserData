#!/usr/bin/env python

import sys

inputData = list()
for line in sys.stdin:
    tokens = line.rstrip().split('\t')
    inputData.append(tokens)
for col in range(len(tokens)):
    outputRow = ""
    delimiter = ""
    for row in range(len(inputData)):
        outputRow = outputRow + delimiter + str(inputData[row][col])
        delimiter = "\t"
    print outputRow
