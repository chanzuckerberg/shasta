#!/usr/bin/python3

import os
import sys

helpMessage = """
This script will convert a fastq file to fasta.

Invoke with two arguments:
- Name of the input fastq file.
- Name of the output fasta file.
"""

if not len(sys.argv)==3:
    print(helpMessage)
    exit(1) 
    
inputFileName = sys.argv[1]
outputFileName = sys.argv[2]

command = 'cat ' + inputFileName + ' | awk \'{if(NR%4==1) {printf(">%s\\n",substr($0,2));} else if(NR%4==2) print;}\' > ' + outputFileName
    
os.system(command)    

