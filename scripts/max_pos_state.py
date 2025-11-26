#!/usr/bin/env python
import numpy as np
import pandas as pd
import argparse as ap
import jpype, sys, os
import jpype.imports
from jpype.types import *
from scipy.stats import mode

strands = pd.read_csv("strands.txt", header=None, index_col=0)

def extract_gene(df, gene):
    df1 = df[df["Ensembl_geneid"] == gene].sort_values("pos(1-based)").reset_index().copy()
    if gene in strands.index and strands[1][gene] == '-':
        df1 = df1[::-1].reset_index(drop=True)

    unique_poss = df1["pos(1-based)"].unique()
    pos_to_index = {}
    for i in range(len(unique_poss)):
        pos_to_index[unique_poss[i]] = i

    df1["new_pos"] = [pos_to_index[item] for item in df1["pos(1-based)"]]
    return df1


def generate_sample3(genedf_data_ordered, target_row, target_pos, flank_len=5, num_samples=10):
    target_index = flank_len
    start = max(0, target_pos - flank_len)
    if start == 0:
        target_index = target_pos
    end = min(len(genedf_data_ordered), target_pos + flank_len + 1)
    sub_data = genedf_data_ordered[start:end]
    samples = []
    for i in range(num_samples):
        rand_inds = np.random.randint(0, [len(item) for item in sub_data])
        sample = [sub_data[i][rand_inds[i]] for i in range(len(rand_inds))]
        sample[target_index] = target_row
        samples.append(sample)
    return samples, target_index


parser = ap.ArgumentParser()
parser.add_argument('-c', '--chrm', help="Chromosome id", nargs='?')
parser.add_argument('-v', '--input_version', help="Version of dbNSFP database used for input. Default 5.0a",
                    type=str, default="5.0a", const="5.0a", nargs='?')
parser.add_argument('-m', '--model_path', help="Path to HMM model file", type=str)
parser.add_argument('-f', '--flank_len', help="Length of flanking window for sampling. Default 3",
                    type=int, default=3, const=3, nargs='?')
parser.add_argument('-o', '--output', help="Output path. Default: max_states directory under the same path as the model file",
                    type=str, default="", const="", nargs='?')
parser.add_argument('-t', '--threshold_file', help="File containing thresholds", type=str, default="", const="", nargs='?')
args = parser.parse_args()

chrm = args.chrm
flank_len = args.flank_len

df = pd.read_csv("dbNSFP%s/datafiles/new/dbnsfp_chr%s.txt" % (args.input_version, chrm), sep='\t', index_col=False)
df_sorted = df.sort_values(["Ensembl_geneid", "pos(1-based)"]).reset_index()

remove = ["MutationTaster_converted_rankscore", "phastCons470way_mammalian_rankscore", "phastCons100way_vertebrate_rankscore"]
sorted_features = [item for item in df.columns[6:] if ((item not in remove) and ("rankscore" in item))]


if args.threshold_file == "":
    tfile = "max_thresholds_%s.txt" % (args.input_version)
else:
    tfile = args.threshold_file
#thresholds = pd.read_csv("max_thresholds_%s.txt" % (args.input_version), index_col=0)['0']
thresholds = pd.read_csv(tfile, index_col=0)['0']
delta = 0.00001

# Load chromhmm
# model_file = "dbNSFP4.8a/models/all_%d_missing/model_%d.txt" % (num_states, num_states)
model_file = args.model_path
if args.output == "":
    output_path = os.path.dirname(model_file) + "/max_states"
else:
    output_path = args.output
os.makedirs(output_path, exist_ok=True)


jpype.startJVM(classpath=["/u/home/l/luke0321/project-ernst/scoreHMM/ChromHMM/ChromHMM.jar"])
c = jpype.JClass("edu.mit.compbio.ChromHMM.ChromHMM")
chromhmm = c(model_file)

sys.stdout = open(output_path + "/max_states_chr%s.bed" % (chrm), "w+")

genes = df["Ensembl_geneid"].unique()

for g in genes:
    df1 = extract_gene(df, g)

    missing_mask = np.asarray((df1[sorted_features] == -1).values, dtype=int)
    binarized = np.asarray(df1[sorted_features].ge(thresholds - delta, axis=1)[sorted_features].values, dtype=int)
    binarized[missing_mask == 1] = 2

    df1_data = binarized
    df1_pos = df1["new_pos"].values
    df1_data_ordered = [df1_data[df1_pos == i] for i in range(df1_pos[-1] + 1)] # ith row of this is all variants at ith position

    old_pos = list(df1["pos(1-based)"])
    refs = list(df1["ref"])
    alts = list(df1["alt"])
    aachanges = list(df1["AA"])

    for i in range(len(df1_data)):
        samples, target_index = generate_sample3(df1_data_ordered, df1_data[i], df1_pos[i], flank_len=flank_len, num_samples=9)
        states = []
        for sample in samples:
            data = np.asarray(sample)
            states.append(chromhmm.getMaxStateAtPos(data.tolist(), target_index))
        max_state = max(set(states), key = states.count)
        sys.stdout.write("chr%s\t%d\t%d\t%s\t%s\tE%d\t%s\n" % (chrm, old_pos[i] - 1, old_pos[i], refs[i], alts[i], max_state, aachanges[i]))
