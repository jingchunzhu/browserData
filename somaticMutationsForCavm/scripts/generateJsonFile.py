#!/usr/bin/env python

import argparse
import json


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("-a", "--assembly", help="Assembly (genome). Default: hg19", default="hg19")
    parser.add_argument("-c", "--cohort", help="Cohort (no default)")
    parser.add_argument("-d", "--description", help="Description, no default")
    parser.add_argument("-l", "--longTitle", help="Long title, no default.")
    parser.add_argument("-n", "--name", help="Name")
    parser.add_argument("-r", "--reference", help="Reference file pathname")
    args = parser.parse_args()

    jsonData = {"name": args.name, "type": "mutationVector", ":assembly": args.assembly,
                "cohort": args.cohort, ":dataSubType": "SNP_small_INDEL",
                "version": "12-13-2013", "shortTitle": args.cohort + " mutations",
                "longTitle": args.longTitle,
                "description": args.description,
                "reference": args.reference} 
    print json.dumps(jsonData, indent=2)


if __name__ == '__main__':
        main()
    

