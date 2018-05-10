#!/usr/bin/python3

import Nanopore2

readId = int(input('Enter a ReadId: '))
fileName = 'Markers-' + str(readId) + '.csv'

a = Nanopore2.Assembler()
a.accessMarkers()
a.writeMarkers(readId=readId, fileName=fileName)

