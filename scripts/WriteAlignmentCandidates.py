#!/usr/bin/python3
import argparse
import shasta
import sys


def main(useReadName, verbose):
    a = shasta.Assembler()
    print("Accessing data ...")
    a.accessMarkers()
    a.accessAlignmentCandidates()

    if verbose:
        a.accessAlignmentData()
        a.accessReadGraph()

    print("Writing to AlignmentCandidates.csv ...")
    a.writeAlignmentCandidates(useReadName, verbose)
    print("Done")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(prog='WriteAlignmentCandidates',
                                     formatter_class=argparse.RawDescriptionHelpFormatter,
                                     epilog="Program description:\n\tWrite all alignment candidates (pairs of ReadId) "
                                            "to a file named AlignmentCandidates.csv")

    parser.add_argument(
        "--useReadName",
        dest="useReadName",
        action="store_true",
        help="If this boolean flag is specified, the full read name"
             " from the input sequence file will be written instead of the ID"
    )

    parser.add_argument(
        "--verbose",
        dest="verbose",
        action="store_true",
        help="If this boolean flag is specified, the result of filtering the"
             " candidate will also be written to the csv"
    )

    args = parser.parse_args()

    main(
        useReadName=args.useReadName,
        verbose=args.verbose
    )