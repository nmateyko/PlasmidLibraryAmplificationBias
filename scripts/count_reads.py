#!/usr/bin/env python

# Takes a list of sequences and counts the occurrences of each unique sequence. Similar sequences are grouped together to account
# for sequencing errors or mutations during plasmid amplification/library prep.

import gzip
import argparse
from collections import Counter
import numpy as np
import random
from umi_tools import UMIClusterer

parser = argparse.ArgumentParser()
parser.add_argument('-i', dest='inFP',	metavar='<inFile>', help='Input file', required=True)
parser.add_argument('-o', dest='outFP', metavar='<outFile>', help='Where to output results', required=True)
parser.add_argument('-d', dest='downsample', type=int, metavar='<downsample n>', help='Number of reads to downsample to', required=True)
parser.add_argument('-l', dest='logFP', metavar='<logFile>', help='Where to output errors/warnings', required=True)

args = parser.parse_args()
clusterer = UMIClusterer(cluster_method="directional")

with gzip.open(args.inFP, 'rt') as inFile:
  seqs = inFile.read().splitlines()
  
seqs = [seq.encode() for seq in seqs] 
print(len(seqs))
seqs = random.sample(seqs, args.downsample)
counts = dict(Counter(seqs))
clustered = clusterer(counts, threshold=2)
cluster_counts = []
for cluster in clustered:
  sum = 0
  for seq in cluster:
    sum += counts[seq]
  cluster_counts.append(sum)

print(np.array(cluster_counts))
with open(args.outFP, 'wb') as f:
  np.save(f, cluster_counts)

