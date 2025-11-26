#!/usr/bin/env python

import numpy as np
import pandas as pd
import argparse as ap
import os, sys, gzip


# https://stackoverflow.com/questions/736043/checking-if-a-string-can-be-converted-to-float-in-python
def is_float(element: any) -> bool:
    #If you expect None to be passed:
    if element is None:
        return False
    try:
        float(element)
        return True
    except ValueError:
        return False

def process_pos(lines):
    x_alts = []
    selected_lines = []
    for line in lines:
        text = line.split('\t')
        alt = text[3]
        aaref = text[4]
        aaalt = text[5]
        if aaref == 'X' or aaalt == 'X':
            x_alts.append(alt)

    #if len(x_alts) >= 1:
    #    return []

    for line in lines:
        text = line.split('\t')
        alt = text[3]
        aaref = text[4]
        aaalt = text[5]
        if alt in x_alts or aaref == aaalt or aaref == '.' or aaalt == '.':
            continue
        else:
            # Select only missense
            selected_lines.append(line)
    return selected_lines


parser = ap.ArgumentParser()
parser.add_argument('-c', '--chrm', help="Chromosome id", nargs='?')
parser.add_argument('-v', '--input_version', help="Version of dbNSFP database used for input. Default 5.0a", 
                    type=str, default="5.0a", const="5.0a", nargs='?')
parser.add_argument('-o', '--output', help="Output path. Default: dbNSFP[version]/datafiles", 
                    type=str, default="", const="", nargs='?')
parser.add_argument('-g', '--gnomad', help="Output the gnomAD allele frequencies", action="store_true")
parser.add_argument('-l', '--clinvar', help="Output the clinvar consequences", action="store_true")
args = parser.parse_args()


chrm = args.chrm
if args.output == "":
    output_path = "dbNSFP%s/datafiles" % (args.input_version)
else:
    output_path = args.output
os.makedirs(output_path, exist_ok=True)
path = "/u/home/l/luke0321/project-ernst/scoreHMM/dbnsfp/dbNSFP%s/dbNSFP%s_variant.chr%s.gz" % (args.input_version, args.input_version, chrm)
output = output_path + "/dbnsfp_chr%s.txt" % chrm


header_inds = []
with gzip.open(path, "rt") as file, open(output, "w+") as output_file:
    line = file.readline()
    header = line.split('\t')
    # Indices of columns containing rankscores
    header_inds = [i for i in range(len(header)) if "rankscore" in header[i]]
    if args.gnomad:
        header_inds.append(header.index("gnomAD4.1_joint_AF"))
        header_inds.append(header.index("gnomAD4.1_joint_AC"))
    if args.clinvar:
        header_inds.append(header.index("clinvar_clnsig"))
        header_inds.append(header.index("clinvar_trait"))
        header_inds.append(header.index("clinvar_review"))
        
    # Indices 0, 1, 2, 3, 13, correspond to Chr, pos, ref, alt, geneid
    output_file.write('\t'.join([header[0], header[1], header[2], header[3], header[13], "AA"]) + '\t')
    output_file.write('\t'.join([header[i] for i in header_inds]) + '\n')

    line = file.readline()
    prev_pos = -1
    pos_lines = []
    while line:
        text = line.split('\t')
        pos = int(text[1])

        if pos == prev_pos:
            pos_lines.append(line)

        else:
            if prev_pos != -1:
                selected_lines = process_pos(pos_lines)
                for sl in selected_lines:
                    sl_text = sl.split('\t')
                    for i in header_inds:
                        if not is_float(sl_text[i]):
                            # If it is a '.' convert to -1
                            if sl_text[i] == '.':
                                sl_text[i] = '-1'
                            # Otherwise print it out
                            else:
                                print(sl_text[i])

                    selected_cols = [sl_text[i] for i in header_inds]
                    output_file.write('\t'.join([sl_text[0], sl_text[1], sl_text[2], sl_text[3], sl_text[13].split(';')[0], "%s-%s" % (sl_text[4], sl_text[5])] + selected_cols) + '\n')

            prev_pos = pos
            pos_lines = [line]

        line = file.readline()
