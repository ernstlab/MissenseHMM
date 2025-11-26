#!/usr/bin/env python
import numpy as np
import pandas as pd
import argparse as ap
import os, sys, gzip

parser = ap.ArgumentParser()
parser.add_argument('-c', '--chrm', help="Chromosome id", nargs='?')
parser.add_argument('-v', '--input_version', help="Version of dbNSFP database used for input. Default 5.0a", 
                    type=str, default="5.0a", const="5.0a", nargs='?')
parser.add_argument('-i', '--input', help="Input path. Default: dbNSFP[version]/datafiles", 
                    type=str, default="", const="", nargs='?')
parser.add_argument('-o', '--output', help="Output path. Default: dbNSFP[version]/datafiles/samples", 
                    type=str, default="", const="", nargs='?')
parser.add_argument('-t', '--threshold_file', help="File containing thresholds", type=str, default="", const="", nargs='?')
args = parser.parse_args()

chrm = args.chrm

if args.input == "":
    input_path = "dbNSFP%s/datafiles" % (args.input_version)
else:
    input_path = args.input
    
if args.output == "":
    output_path = "dbNSFP%s/datafiles/samples" % (args.input_version)
else:
    output_path = args.output
    
if args.threshold_file == "":
    tfile = "max_thresholds_%s.txt" % (args.input_version)
else:
    tfile = args.threshold_file
    
os.makedirs(output_path, exist_ok=True)


df = pd.read_csv("%s/dbnsfp_chr%s.txt" % (input_path, chrm), sep='\t', index_col=False)
thresholds = pd.read_csv(tfile, index_col=0)['0']
delta = 0.00001

df = df.groupby('pos(1-based)').apply(lambda x: x.sample(1)).droplevel(1)
df["start"] = df.index - 1
df["end"] = df.index
df["chr"] = "chr%s" % chrm

df = df[["chr", "start", "end", "ref", "alt", "Ensembl_geneid"] + list(df.columns[6:-3])]

df_sorted = df.sort_values(["Ensembl_geneid", "start"])
genes = list(df_sorted["Ensembl_geneid"])
new_starts = np.zeros(len(df_sorted), dtype=int)
current_cnt = 0
prev_gene = genes[0]
for i in range(1, len(genes)):
    current_gene = genes[i]
    if current_gene != prev_gene:
        current_cnt = 0
        new_starts[i] = current_cnt
    else:
        current_cnt += 1
        new_starts[i] = current_cnt
    prev_gene = current_gene

df_sorted["start_g"] = new_starts
df_sorted["end_g"] = new_starts + 1

test = df_sorted.head().copy()

remove = ["MutationTaster_converted_rankscore", "phastCons470way_mammalian_rankscore", "phastCons100way_vertebrate_rankscore"]
features = [item for item in df.columns[6:] if ((item not in remove) and ("rankscore" in item))]

strands = pd.read_csv("strands.txt", header=None, index_col=0)

for name, group in df_sorted.groupby("Ensembl_geneid"):
    missing_mask = np.asarray((group[features] == -1).values, dtype=int)
    data = np.asarray(group[features].ge(thresholds - delta, axis=1)[features].values, dtype=int)
    data[missing_mask == 1] = 2
    
    strand = strands[1][name]
    if strand == '-':
        data = np.flip(data, axis=0)
    
    with open(output_path + "/%s_binary.txt" % name, "w+") as file:
        file.write("cell1\t%s\n" % name)
        file.write('\t'.join(features) + '\n')
        for line in data:
            file.write('\t'.join([str(item) for item in line]) + '\n')