#!/usr/bin/env python
import numpy as np
import sys, glob

score_index = int(sys.argv[1])
# chrms = [str(item) for item in range(1, 23)] + ['X']
threshold = 90

input_dir = sys.argv[2]
output_dir = sys.argv[3]
files = glob.glob(input_dir + "/dbnsfp_chr*.txt")
total_variants = 0
for f in files:
    with open(f) as file:
        total_variants += (sum(1 for _ in file) - 1)
        
scores = np.full(total_variants, -1, dtype=float)
cnt = 0
score_name = ""

for f in files:
    with open(f) as file:
        line = file.readline()
        text = line.split()
        score_name = text[6 + score_index]
        line = file.readline()
        while line:
            text = line.split()
            scores[cnt] = float(text[6 + score_index])
            cnt += 1
            line = file.readline()

with open("%s/dbnsfp_%d.log" % (output_dir, score_index), "w+") as file:
    file.write("score\tnonmissing\t90p\tones\n")
    nonmissing = scores[scores != -1]
    sp = np.percentile(nonmissing, threshold)
    file.write("%s\t%d\t%f\t%d\n" % (score_name, len(nonmissing), sp, len(nonmissing[nonmissing >= sp])))