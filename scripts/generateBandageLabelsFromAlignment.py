from collections import defaultdict
import argparse


def getChromosomeColor(refName):
    if refName == "chr1":
        return "sea green"
    elif refName == "chr2":
        return "tomato"
    elif refName == "chr3":
        return "dark orange"
    elif refName == "chr4":
        return "dark orchid"
    elif refName == "chr5":
        return "firebrick"
    elif refName == "chr6":
        return "deep pink"
    elif refName == "chr7":
        return "#800000"
    elif refName == "chr8":
        return "web green"
    elif refName == "chr9":
        return "green yellow"
    elif refName == "chr10":
        return "gold"
    elif refName == "chr11":
        return "midnight blue"
    elif refName == "chr12":
        return "light coral"
    elif refName == "chr13":
        return "light green"
    elif refName == "chr14":
        return "red"
    elif refName == "chr15":
        return "aquamarine"
    elif refName == "chr16":
        return "dodger blue"
    elif refName == "chr17":
        return "orange red"
    elif refName == "chr18":
        return "yellow"
    elif refName == "chr19":
        return "teal"
    elif refName == "chr20":
        return "aqua"
    elif refName == "chr21":
        return "hot pink"
    elif refName == "chr22":
        return "light salmon"
    elif refName == "chrX":
        return "yellow green"
    elif refName == "chrY":
        return "olive drab"
    else:
        print("Warning using gray for unknown chromosome: " + refName)
        return "gray"


def main(path):
    primaryAlignments = defaultdict(lambda: defaultdict(int))

    with open(path, 'r') as file:
        for l,line in enumerate(file):
            tokens = line.strip().split('\t')

            queryName = tokens[0]
            refName = tokens[5]
            length = int(tokens[10])
            
            if not refName.startswith("chr") or len(refName) > 5:
                continue

            # Store all alignment lengths in a dictionary
            primaryAlignments[queryName][refName] += length

    with open("GfaLabels.csv", 'w') as file:
        file.write("name,color,refName\n")

        for queryName,alignments in primaryAlignments.items():
            maxLength = 0
            bestRef = None

            # Find the reference contig with the most aligned length
            for refName,length in alignments.items():
                print(queryName, refName, length)
                
                if length > maxLength:
                    maxLength = length
                    bestRef = refName
                    
            file.write(queryName)
            file.write(',')
            file.write(getChromosomeColor(bestRef))
            file.write(',')
            file.write(bestRef)
            file.write('\n')


if __name__ == "__main__":
    parser = argparse.ArgumentParser(prog='generateBandageLabelsFromAlignment',
                                     formatter_class=argparse.RawDescriptionHelpFormatter,
                                     epilog="Program description:\n\tuse a PAF alignment tsv to generate a  "
                                            "karyotype-style coloring and labeling for Bandage")

    parser.add_argument(
        "--input","-i",
        required=True,
        type=str,
        help="Path to any PAF formatted TSV that was generated from minimap2/winnowmap for the shasta contigs."
             "Only chromosomes with names chr1-22 or chrX/Y will be colored"
    )

    args = parser.parse_args()

    main(
        path=args.input
    )
