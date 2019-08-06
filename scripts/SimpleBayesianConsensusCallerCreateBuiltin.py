#!/usr/bin/python3

import argparse
import os


def main(inputPath, outputPath):

    is_header = False
    with open(inputPath, "r") as input_file, open(outputPath, "w") as output_file:
        for line in input_file:
            line = line.strip()
            comma = ""

            if line.startswith(">"):
                header = line[1:]
                header_tokens = header.split(" ")
                is_header = True

            elif len(line) > 0:
                if header_tokens[-1] == "Name":
                    if is_header:
                        output_file.write("configurationName = ")
                        is_header = False

                    output_file.write("\"" + line + "\";\n")

                elif header_tokens[-1] == "prior":
                    if is_header and header_tokens[0] == "AT":
                        output_file.write("priors = {{")
                        is_header = False
                    else:
                        comma = ","

                    output_file.write(comma+"\n    {" + line + "}")

                elif header_tokens[-1] == "likelihood":
                    if is_header and header_tokens[0] == "A":
                        output_file.write("\n}};\n")
                        output_file.write("probabilityMatrices = {{\n{")
                        is_header = False
                    elif is_header and header_tokens[0] != "A":
                        output_file.write("\n},\n{")
                        is_header = False
                    else:
                        comma = ","

                    output_file.write(comma+"\n    {" + line + "}")

        output_file.write("\n}}};\n")


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--input",
        type=str,
        required=True,
        help="Path to SimpleBayesianConsensusCaller configuration file."
    )
    parser.add_argument(
        "--output",
        type=str,
        required=True,
        help="Path to .hpp file to be created."
    )

    args = parser.parse_args()

    main(inputPath=args.input, outputPath=args.output)
