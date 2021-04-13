#!/usr/bin/python3
import argparse
import shasta
import sys


def main(useReadName):
    a = shasta.Assembler()
    print("Accessing data ...")
    a.accessMarkers()
    a.accessAlignmentCandidates()
    print("Writing to AlignmentCandidates.csv ...")
    a.writeAlignmentCandidates(useReadName)
    print("Done")


if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    parser.add_argument(
        "--useReadName",
        dest="useReadName",
        action="store_true",
        help="If this boolean flag is specified, the full read name"
             " from the input sequence file will be written instead of the ID"
    )

    args = parser.parse_args()

    main(
        useReadName=args.useReadName
    )