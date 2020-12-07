#!/usr/bin/env bash
set -e

./make_test_fastas.py
fastaq to_perfect_reads ref_for_sim_reads.fa reads.fq.gz 200 1 50 50

bwa index ref.fa
bwa mem ref.fa reads.fq.gz | samtools sort -O bam -o mapped_reads.bam
rm ref_for_sim_reads.fa
fastaq deinterleave reads.fq.gz reads_1.fq.gz reads_2.fq.gz
rm reads.fq.gz
