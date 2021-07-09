#!/usr/bin/python3
from collections import defaultdict
from sklearn import metrics
import argparse
import os


class PhasedRead:
    def __init__(self, id, name, component, phase):
        self.id = id
        self.name = name
        self.component = int(component)
        self.phase = int(phase)


# Iterate the shasta CSV and yield one PhasedRead data object at a time
def parseShastaPhaseCsv(shastaPath):
    with open(shastaPath, 'r') as file:
        for l,line in enumerate(file):
            if l == 0:
                # Skip header line
                continue

            # OrientedReadId,ReadName,Component,Phase
            tokens = line.strip().split(',')

            # Ignore unphased reads (for now)
            if (len(tokens[-1]) == 0) and (len(tokens[-2]) == 0):
                continue

            yield PhasedRead(tokens[0], tokens[1], tokens[2], tokens[3])


def main(hap1MarginPath, hap2MarginPath, shastaPath):
    discordantCsvPath = os.path.join(os.path.dirname(shastaPath), "DiscordantPhasedReads.csv")

    marginPhases = (set(),set())
    shastaComponents = defaultdict(lambda: (set(),set()))
    nameToId = dict()

    # Read all the Shasta phase assignments
    for read in parseShastaPhaseCsv(shastaPath):
        shastaComponents[read.component][read.phase].add(read.name)
        nameToId[read.name] = read.id

        if read.name in shastaComponents[read.component][1-read.phase]:
            exit("ERROR: read contained both in hap1 and hap2 of Shasta phases: " + read.name)

    # Read the Margin phase assignments
    with open(hap1MarginPath, 'r') as margin0File, open(hap2MarginPath, 'r') as margin1File:
        for line in margin0File:
            name = line.strip()
            marginPhases[0].add(name)

        for line in margin1File:
            name = line.strip()
            marginPhases[1].add(name)

    # Compare Shasta to Margin and write results to CSV + stdout
    with open(discordantCsvPath, 'w') as file:
        file.write("OrientedReadId,ReadName,Component,ShastaPhase,MarginPhase\n")

        for componentId,component in shastaComponents.items():
            print("Evaluating component:", componentId)

            names = list()
            shastaPhaseArray = list()
            marginPhaseArray = list()

            # Compare phase assignments for Shasta and the user-provided dataset, for each read in each component
            for p,phase in enumerate(component):
                for readName in phase:
                    marginPhase = -1

                    if readName in marginPhases[0]:
                        marginPhase = 0

                    if readName in marginPhases[1]:
                        if marginPhase == 0:
                            exit("ERROR: read contained both in hap1 and hap2 of Margin phases: " + readName)
                        else:
                            marginPhase = 1

                    # Only store results for reads that also exist in the user-provided dataset
                    if marginPhase >= 0:
                        names.append(readName)
                        shastaPhaseArray.append(p)
                        marginPhaseArray.append(marginPhase)

            countMatrix = [[0,0],[0,0]]
            confusionMatrix = [[0,0],[0,0]]

            for i in range(len(names)):
                countMatrix[0][shastaPhaseArray[i]] += 1
                countMatrix[1][marginPhaseArray[i]] += 1

                confusionMatrix[shastaPhaseArray[i]][marginPhaseArray[i]] += 1

            # Compute the weight of the diagonal and antidiagonal, and infer which is True Positives
            diagonalSum = confusionMatrix[0][0] + confusionMatrix[1][1] + 1e-12
            antidiagonalSum = confusionMatrix[1][0] + confusionMatrix[0][1]

            if 2 > (antidiagonalSum / diagonalSum) > 0.5:
                print("WARNING: ratio of antidiagonalSum:diagonalSum is only " + str(antidiagonalSum / diagonalSum))

            # Indicate which coordinates in the confusion matrix represent discordant reads between shasta's
            # and the user-provided phase assignments
            discoordinateA = (0,1)
            discoordinateB = (1,0)

            # Iterate reads again and write the ones that are discordant with the user-provided phases
            if antidiagonalSum > diagonalSum:
                discoordinateA = (0,0)
                discoordinateB = (1,1)

            for i in range(len(names)):
                coordinate = (shastaPhaseArray[i], marginPhaseArray[i])

                if coordinate == discoordinateA or coordinate == discoordinateB:
                    file.write(nameToId[names[i]])
                    file.write(',')
                    file.write(names[i])
                    file.write(',')
                    file.write(str(componentId))
                    file.write(',')
                    file.write(str(shastaPhaseArray[i]))
                    file.write(',')
                    file.write(str(marginPhaseArray[i]))
                    file.write('\n')

            print("Shasta phase counts:",countMatrix[0])
            print("Margin phase counts:",countMatrix[1])

            print("Confusion matrix:")
            print(confusionMatrix[0])
            print(confusionMatrix[1])

            adjustedRandIndex = metrics.adjusted_rand_score(marginPhaseArray, shastaPhaseArray)
            print("ARI:", adjustedRandIndex)
            print()


if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    parser.add_argument(
        "-1","--hap1",
        type=str,
        required=True,
        help="Path of Marginphase txt file for hap1, containing one read name per line"
    )
    parser.add_argument(
        "-2","--hap2",
        type=str,
        required=True,
        help="Path of Marginphase txt file for hap2, containing one read name per line"
    )
    parser.add_argument(
        "-s","--shasta",
        type=str,
        required=True,
        help="Path of Shasta phase components csv (Phasing.csv)"
    )

    args = parser.parse_args()

    main(
        hap1MarginPath=args.hap1,
        hap2MarginPath=args.hap2,
        shastaPath=args.shasta,
    )

