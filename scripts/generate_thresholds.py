#!/usr/bin/env python
import numpy as np
import pandas as pd
import glob, gzip, os, sys


logs = glob.glob(sys.argv[1] + "/*.log")
output = sys.argv[2]

data = []
for l in logs:
    data.append(pd.read_csv(l, sep='\t'))

df = pd.concat(data, axis=0).reset_index(drop=True)
df.index = df["score"]

df["totals"] = 70842617 - 23
df["missing_rate"] = 1 - (df["nonmissing"] / df["totals"])
df["ones_rate"] = df["ones"] / df["nonmissing"]
df[["90p", "missing_rate", "ones_rate"]].sort_values("ones_rate", ascending=False)[:10]
df[["90p", "missing_rate", "ones_rate"]].sort_values("missing_rate", ascending=False)[:10]
output_df = pd.DataFrame(df["90p"].copy())
output_df.columns = [0]
output_df.to_csv(output + "/thresholds.txt")