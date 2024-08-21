#!/usr/bin/env python3

# This script will produce summaries of various NLR types from NLR Annotator

# Import python modules

import argparse

# Prepare function to parse command line arguments


def parse_args():
    parser = argparse.ArgumentParser(description='Summarise NLR Annotator \
        results')
    parser.add_argument('--input', required=True, type=str, help='Path to \
        NLR Annotator text file to summarise')
    parser.add_argument('--output', required=True, type=str, help='Path to \
        summary file')

    return parser.parse_args()

# Prepare function to count different gene types


def load_file(input_file: str):
    lines = open(input_file).readlines()
    contig_set = set()
    count = 0
    pseudogenes = 0
    genes = 0
    complete = 0
    complete_pseudogenes = 0
    for line in lines:
        count += 1
        line = line.rstrip()
        split_line = line.split('\t')
        nlr_type = split_line[2]
        contig = split_line[0]
        contig_set.add(contig)
        if nlr_type == "complete (psuedogene)" or nlr_type == "partial (pseudogene)":
            pseudogenes += 1
        if nlr_type == "complete" or nlr_type == "partial":
            genes += 1
        if nlr_type == "complete":
            complete += 1
        if nlr_type == "complete (psuedogene)":
            complete_pseudogenes += 1
    contig_count = len(contig_set)

    return contig_count, count, pseudogenes, genes, complete, complete_pseudogenes

# Prepare function to write out summary file


def write_output(input_file: str, output_file: str):
    contig_count, count, pseudogenes, genes, complete, complete_pseudogenes = load_file(input_file)
    with open(output_file, "w") as out_file:
        header_list = ["NLR Contigs", "NLR Count", "Pseudogenous NLRs", "NLR Genes", "Complete NLRs", "Complete Pseudogenous NLRs"]
        header_string = "\t".join(header_list)
        out_file.write(header_string)
        out_file.write("\n")
        list_to_write = [str(contig_count), str(count), str(pseudogenes), str(genes), str(complete), str(complete_pseudogenes)]
        string_to_write = "\t".join(list_to_write)
        out_file.write(string_to_write)
        out_file.write("\n")

# Prepare main function


def main():
    args = parse_args()
    input_file = args.input
    output_file = args.output
    write_output(input_file, output_file)
    

if __name__ == '__main__':
    main()
