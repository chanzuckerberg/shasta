#!/usr/bin/python3

import shasta
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('readName', type=str)
arguments = parser.parse_args()

a = shasta.Assembler()
print(a.getReads().getReadId(arguments.readName))

