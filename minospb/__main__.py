#!/usr/bin/env python3

import argparse
import logging

import minospb


def main(args=None):
    parser = argparse.ArgumentParser(
        prog="minospb",
        usage="minospb <command> <options>",
        description="minospb: variant call adjudication",
    )

    parser.add_argument("--version", action="version", version=minospb.__version__)

    subparsers = parser.add_subparsers(title="Available commands", help="", metavar="")

    # -------------------------- make_mask ------------------------------------
    subparser_make_mask = subparsers.add_parser(
        "make_mask",
        help="Make genome mask, where low read depth or read consensus not good enough",
        usage="minospb make_mask [options] <ref_fa> <reads1> <reads2> <outdir>",
        description="Make genome mask, where low read depth or read consensus not good enough",
    )
    subparser_make_mask.add_argument(
        "--min_depth",
        type=int,
        help="Minimum read depth [%(default)s]",
        metavar="INT",
        default=5,
    )
    subparser_make_mask.add_argument(
        "--min_depth_pc",
        type=float,
        help="Minimum percent of reads that must agree with reference [%(default)s]",
        metavar="FLOAT",
        default=90,
    )
    subparser_make_mask.add_argument("ref_fa", help="FASTA file of reference genome")
    subparser_make_mask.add_argument("reads1", help="FASTQ file of forwards reads")
    subparser_make_mask.add_argument("reads2", help="FASTQ file of reverse reads")
    subparser_make_mask.add_argument(
        "outdir", help="Output directory (will be created)"
    )
    subparser_make_mask.set_defaults(func=minospb.tasks.make_mask.run)

    # ---------------------------- index_ref_genome ---------------------------
    subparser_index_ref_genome = subparsers.add_parser(
        "index_ref_genome",
        help="Make index files for a reference genome",
        usage="minospb index_ref_genome <fasta_in> <outdir>",
        description="Make index files for a reference genome",
    )
    subparser_index_ref_genome.add_argument(
        "fasta_in", help="Name of FASTA file to be indexed",
    )
    subparser_index_ref_genome.add_argument(
        "outdir", help="Name of output directory (must not already exist)",
    )
    subparser_index_ref_genome.set_defaults(func=minospb.tasks.index_ref_genome.run)

    # ---------------------------- run_one_sample -----------------------------
    subparser_run_one_sample = subparsers.add_parser(
        "run_one_sample",
        help="Variant call and run all benchmarking tools on one sample",
        usage="minospb run_one_sample [options] <sample_name> <sample_dir> <ref_dir> <truth_fasta> <reads1> <reads2>",
        description="Variant call and run all benchmarking tools on one sample",
    )
    subparser_run_one_sample.add_argument(
        "--clockwork",
        action="store_true",
        help="Just run the 'clockwork' pipeline, which means trim reads, map, call with cortex+samtools, and then run minos",
    )
    subparser_run_one_sample.add_argument(
        "--minos_ref_splits",
        type=int,
        help="Number of reference split chunks when running minos (is passed to minos adjudicate --total_splits) [%(default)s]",
        default=1,
        metavar="INT",
    )
    subparser_run_one_sample.add_argument(
        "--varifier_no_maxmatch",
        action="store_true",
        help="Use the --no_maxmatch option with varifier",
    )
    subparser_run_one_sample.add_argument(
        "--truth_vcf", help="VCF file of truth variant calls", metavar="FILENAME",
    )
    subparser_run_one_sample.add_argument(
        "--ref_mask_bed",
        help="BED file of regions to mask from ref genome when evaluating accuracy of calls",
        metavar="FILENAME",
    )
    subparser_run_one_sample.add_argument(
        "--truth_mask_bed",
        help="BED file of regions to mask from truth genome when evaluating accuracy of calls",
        metavar="FILENAME",
    )
    subparser_run_one_sample.add_argument(
        "--ram",
        type=float,
        help="RAM limit (for tools that have the option) in GB [%(default)s]",
        default=7,
        metavar="FLOAT",
    )
    subparser_run_one_sample.add_argument(
        "sample_name", help="Name of sample",
    )
    subparser_run_one_sample.add_argument(
        "sample_dir", help="Root directory of sample",
    )
    subparser_run_one_sample.add_argument(
        "ref_dir",
        help="Directory with reference files. Made by `minospb index_ref_genome` or by `clockwork reference_prepare`",
    )
    subparser_run_one_sample.add_argument(
        "truth_fasta", help="FASTA file of truth genome",
    )
    subparser_run_one_sample.add_argument(
        "reads1", help="FASTQ file of forwards reads",
    )
    subparser_run_one_sample.add_argument(
        "reads2", help="FASTQ file of reverse reads",
    )
    subparser_run_one_sample.set_defaults(func=minospb.tasks.run_one_sample.run)

    args = parser.parse_args()

    logging.basicConfig(
        format=f"[%(asctime)s minospb %(levelname)s] %(message)s",
        datefmt="%Y-%m-%dT%H:%M:%S",
    )
    log = logging.getLogger()
    log.setLevel(logging.INFO)

    if hasattr(args, "func"):
        args.func(args)
    else:
        parser.print_help()


if __name__ == "__main__":
    main()
