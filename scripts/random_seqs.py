import random
import gzip

LEN = 20
NUM = 10

bases = ['A', 'C', 'G', 'T']
unique_seqs = ["".join(random.choices(bases, k=LEN)) for i in range(NUM)]
unique_counts = range(1, NUM + 1)
mutated_seqs = []

for i, seq in zip(unique_counts, unique_seqs):
  for j in range(i):
    mut_idx = random.choice(range(LEN))
    mut_seq = list(seq)
    mut_base = mut_seq[mut_idx]
    mut_seq[mut_idx] = random.choice(list(set(bases) - set(mut_base)))
    mutated_seqs.append("".join(mut_seq) + '\n')

unique_seqs = [seq + '\n' for seq in unique_seqs]
mutated_seqs = mutated_seqs + unique_seqs

with gzip.open('test_seqs.txt', 'wt') as f:
  f.writelines(mutated_seqs)

with open('true_counts.txt', 'w') as f:
  for i, seq in zip(unique_counts, unique_seqs):
    f.write(f'{seq} {i}\n')