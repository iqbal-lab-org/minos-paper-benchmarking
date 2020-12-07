import os
import subprocess

from clockwork import read_map, read_trim
import pyfastaq
import pysam

from minospb import utils


def pysam_pileup_list_is_ok(pileup_list, ref_base, min_depth, min_depth_percent):
    if len(pileup_list) == 0:
        return False

    counts = {x: 0 for x in ["A", "C", "G", "T", "I", "D"]}
    for s in pileup_list:
        if "+" in s:
            counts["I"] += 1
        elif "-" in s or s == "*":
            counts["D"] += 1
        else:
            counts[s.upper()] += 1

    try:
        ref_cov = counts[ref_base.upper()]
    except KeyError:
        return False

    if ref_cov < min_depth:
        return False

    return 100 * ref_cov / sum(counts.values()) >= min_depth_percent


def bam_to_mask(bam, ref_fasta, outfile, min_depth, min_depth_percent):
    # Reminder of bed coords:
    # Start is 0-based. End is not included.
    # ie the same as python slices.
    # eg 0-10 means mask the first 10 bases
    mask = {}  # ref name -> list of start/end coords
    current_ref = None
    start = None
    end = None
    ref_seqs = {}
    pyfastaq.tasks.file_to_dict(ref_fasta, ref_seqs)
    samfile = pysam.AlignmentFile(bam, "rb")

    with open(outfile, "w") as f_out:
        for iterator in samfile.pileup():
            pileup = iterator.get_query_sequences(add_indels=True)
            ref_name = iterator.reference_name
            ref_pos = iterator.reference_pos
            ref_base = ref_seqs[ref_name][ref_pos]
            pileup_ok = pysam_pileup_list_is_ok(pileup, ref_base, min_depth, min_depth_percent)

            if pileup_ok:
                if current_ref is not None:
                    print(current_ref, start, end + 1, sep="\t", file=f_out)
                    current_ref = None
                    start = None
                    end = None
                continue

            if current_ref == ref_name and ref_pos == end + 1:
                end += 1
            else:
                if current_ref is not None:
                    print(current_ref, start, end + 1, sep="\t", file=f_out)
                current_ref = ref_name
                start = ref_pos
                end = ref_pos

        if current_ref is not None:
            print(current_ref, start, end + 1, sep="\t", file=f_out)


def trim_and_map_reads_and_make_mask(
    ref_fa, reads1, reads2, outdir, min_depth, min_depth_percent
):
    outdir = os.path.abspath(outdir)
    ref_fa = os.path.abspath(ref_fa)
    reads1 = os.path.abspath(reads1)
    reads2 = os.path.abspath(reads2)
    trim_reads1 = "trimmed_reads_1.fq.gz"
    trim_reads2 = "trimmed_reads_2.fq.gz"
    rmdup_bam = "rmdup.bam"
    mask_bed = "mask.bed"

    original_dir = os.getcwd()
    os.mkdir(outdir)
    os.chdir(outdir)

    read_trim.run_trimmomatic(
        reads1,
        reads2,
        trim_reads1,
        trim_reads2,
        qual_trim="LEADING:10 TRAILING:10 SLIDINGWINDOW:4:15",
    )

    read_map.map_reads(
        ref_fa,
        os.path.abspath(trim_reads1),
        os.path.abspath(trim_reads2),
        rmdup_bam,
        rmdup=True,
    )
    utils.syscall(f"samtools index {rmdup_bam}")

    bam_to_mask(rmdup_bam, ref_fa, mask_bed, min_depth, min_depth_percent)
    os.chdir(original_dir)
    return os.path.join(outdir, mask_bed)
