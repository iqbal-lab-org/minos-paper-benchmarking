#!/usr/bin/env python3

import copy
import random
import pyfastaq
import subprocess

random.seed(42)

seqs = [[random.choice(["A", "C", "G", "T"]) for _ in range(600)]]
seqs.append(copy.copy(seqs[-1]))
seqs.append(copy.copy(seqs[-1]))
seqs[0][200] = "A"
seqs[1][200] = "A"
seqs[2][200] = "G"
seqs[0][201] = "A"
seqs[1][201] = "A"
seqs[2][201] = "G"
seqs[0][400] = "A"
seqs[1][400] = "G"
seqs[2][400] = "G"


sim_reads_fa = "ref_for_sim_reads.fa"
ref_fa = "ref.fa"
with open(sim_reads_fa, "w") as f_sim, open(ref_fa, "w") as f_ref:
    for i, seq in enumerate(seqs):
        fa = pyfastaq.sequences.Fasta(f"seq.{i+1}", "".join(seq))
        print(fa, file=f_sim)
        if i == 0:
            print(fa, file=f_ref)

