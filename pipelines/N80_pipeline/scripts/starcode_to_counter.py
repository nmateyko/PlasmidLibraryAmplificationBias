import csv
import pickle
from collections import Counter
from itertools import islice

# taken from https://gist.github.com/jakebiesinger/759018/1b7d6bd6967780a8bbae743760c37885bdf86467
def read_fastq(fastqfile):
    '''Parse a fastq-formatted file, yielding a (header, sequence, quality) tuple'''
    fastqiter = (l.strip('\n') for l in fastqfile)  # strip trailing newlines
    fastqiter = filter(lambda l: l, fastqiter)  # skip blank lines
    while True:
        fqlines = list(islice(fastqiter, 4))
        if len(fqlines) == 4:
            header1, seq, header2, qual = fqlines
        elif len(fqlines) == 0:
            return
        else:
            raise EOFError("Failed to parse four lines from fastq file!")

        if header1.startswith('@') and header2.startswith('+'):
            yield header1, seq, qual
        else:
            raise ValueError("Invalid header lines: %s and %s for seq %s" % (header1, header2, seq))

def starcode_clusters_to_counter(clustered_path, unclustered_path):
    with open(unclustered_path, 'r') as f:
        reads = [seq for _, seq, _ in read_fastq(f)]
    with open(clustered_path, 'r') as f:
        clustered_lines = [tuple(line) for line in csv.reader(f, delimiter='\t')]

    cluster_counts = {}

    for cluster in clustered_lines:
        indices = [int(i) for i in cluster[3].split(',')]
        cluster_seqs = [reads[i - 1] for i in indices] # Sequence IDs start at 1
        cluster_counts[cluster[0]] = Counter(cluster_seqs)

    return cluster_counts

fastq_path = snakemake.input[0]
cluster_path = snakemake.input[1]
cluster_counter = starcode_clusters_to_counter(cluster_path, fastq_path)

with open(snakemake.output[0], 'wb') as f:
    pickle.dump(cluster_counter, f, protocol=pickle.HIGHEST_PROTOCOL)