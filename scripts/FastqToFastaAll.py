#!/usr/bin/python3

import glob
import os
import sys

helpMessage = """
This script will convert all fastq files in the current directory to fasta.
Invoke without arguments.
"""

if not len(sys.argv)==1:
    print(helpMessage)
    exit(1)
    

# Loop over fastq files in the current directory.  
for fastqFileName in glob.glob('*.fastq'):
    print(fastqFileName)
    fastaFileName = fastqFileName[:-1] + 'a'
    command = 'cat ' + fastqFileName + ' | awk \'{if(NR%4==1) {printf(">%s\\n",substr($0,2));} else if(NR%4==2) print;}\' > ' + fastaFileName   
    os.system(command)        

