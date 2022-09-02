#!/usr/bin/env python

# Takes a list of sequences and counts the occurrences of each unique sequence. Similar sequences are grouped together to account
# for sequencing errors or mutations during plasmid amplification/library prep.

import gzip
import sys
import argparse
from collections import Counter
import matplotlib.pyplot as plt
from umi_tools import UMIClusterer

parser = argparse.ArgumentParser()
parser.add_argument('-i', dest='inFP',	metavar='<inFile>', help='Input file', required=True)
parser.add_argument('-o', dest='outFP', metavar='<outFile>', help='Where to output results', required=True)
parser.add_argument('-l', dest='logFP', metavar='<logFile>', help='Where to output errors/warnings', required=True)
parser.add_argument('-v', dest='verbose', action='count', help='Verbose output?', required=False, default=0)

args = parser.parse_args()
clusterer = UMIClusterer(cluster_method="directional")

with gzip.open(args.inFP, 'rt') as inFile:
  seqs = inFile.read().splitlines() 
  print(len(seqs))
  counts = dict(Counter(seqs))
  print(len(counts))
  clustered = clusterer(counts, threshold=2)
  print(len(clustered))