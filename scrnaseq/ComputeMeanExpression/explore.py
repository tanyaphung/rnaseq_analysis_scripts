import sys
import os
import pandas as pd
import numpy as np
import anndata
from scipy.sparse import csc_matrix
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('-f', '--adata', required=True)
parser.add_argument('-o', '--out', required=True)
args = parser.parse_args()

adata = anndata.read(args.adata)

outfile = open(args.out, "w")

print("### Viewing the columns of the obs.", file=outfile)
print(adata.obs.columns, file=outfile)

print("### View the unique values.", file=outfile)
for i in adata.obs.columns:
    print(f"Find unique values for column: {i}", file=outfile)
    for j in adata.obs[i].unique():
        print(j, file=outfile)

print("### View the adata.var:", file=outfile)
print(adata.var.head(), file=outfile)
print("Size of adata.var:", file=outfile)
print(adata.var.shape, file=outfile)

print("### Check which layer stores the raw count:", file=outfile)
try:
    print("Trying to print adata.X.A[1:25, 1:25]", file=outfile)
    print(adata.X.A[1:25, 1:25], file=outfile)
except: 
    print("command adata.X.A[1:25, 1:25] failed", file=outfile)

try:
    print("Trying to print adata.X[1:25, 1:25]", file=outfile)
    print(adata.X[1:25, 1:25], file=outfile)
except: 
    print("command adata.X[1:25, 1:25] failed", file=outfile)

try:
    print("Trying to print adata.raw.X.A[1:25, 1:25]", file=outfile)
    print(adata.raw.X.A[1:25, 1:25], file=outfile)
except: 
    print("command adata.raw.X.A[1:25, 1:25] failed", file=outfile)

print("### View the adata.obs:", file=outfile)
print(adata.obs.head(), file=outfile)
print("Size of adata.obs:", file=outfile)
print(adata.obs.shape, file=outfile)

print("### Number of female donors: ", file=outfile)
print(len(adata.obs[adata.obs["sex"]=="female"]["donor_id"].unique()), file=outfile)
print("### Number of male donors: ", file=outfile)
print(len(adata.obs[adata.obs["sex"]=="male"]["donor_id"].unique()), file=outfile)