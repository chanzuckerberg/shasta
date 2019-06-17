#!/usr/bin/python3

import argparse
import glob
import json
import collections


helpMessage = """
Creates a single csv file containing summaries for the Shasta assemblies
present in the current directory.
"""

parser = argparse.ArgumentParser(description = helpMessage)
arguments = parser.parse_args()

# Find the json files containing summaries for the individual assemblies.
jsonFileNames = glob.glob('*/AssemblySummary.json')
if not jsonFileNames:
    print('No Shasta assemblies found in the current directory.\n'
        'Invoke this script from a directory containing one or more\n'
        'Shasta assemblies.')
    exit(1)

# Gather the jsons for all the assemblies.
jsons = []
for jsonFileName in jsonFileNames:    
    jsons.append(json.load(open(jsonFileName, 'r'), object_pairs_hook=collections.OrderedDict))

# Open the output file and write the header line.
out = open('AssembliesSummary.csv', 'w')
out.write('Assembly,',)
for jsonFileName in jsonFileNames:
    assemblyName = jsonFileName[:jsonFileName.find('/')]
    out.write('%s,' % assemblyName)
out.write('\n')



# Loop over sections.
firstAssemblyJson = jsons[0]
for sectionName in firstAssemblyJson:
    if sectionName == 'Comment':
        continue
    firstAssemblySection = firstAssemblyJson[sectionName]
    
    # Special treatment of section "Reads discarded on input"
    # which contains reads and bases for each item.
    if sectionName == 'Reads discarded on input':
        for itemName in firstAssemblySection:
           for x in json[sectionName][itemName]:
                out.write('"%s: %s: %s",' % (sectionName, itemName, x))
                for json in jsons:
                    out.write('%s,' % json[sectionName][itemName][x])
                out.write('\n')
        continue
    
    # Loop over items in this section.
    for itemName in firstAssemblySection:
        out.write('"%s: %s",' % (sectionName, itemName))
        
        # Write the values for all the assemblies.
        for json in jsons:
            out.write('%s,' % json[sectionName][itemName])
        
        out.write('\n')
    

