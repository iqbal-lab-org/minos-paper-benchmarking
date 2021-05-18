import glob
import json
import os

import pyfastaq
from cluster_vcf_records import vcf_record

from minospb import utils


def _fix_sv_call(record, ref_seqs):
    assert record.REF == "N"
    if record.ALT[0].startswith("<DEL:"):
        # deleted sequence isn't in the record, so get it from
        #  the reference sequence
        del_length = int(record.INFO["SVLEN"])
        record.REF = ref_seqs[record.CHROM][
            record.POS : record.POS + del_length + 1
        ]
        record.ALT = [record.REF[0]]
    else:
        assert record.ALT[0].startswith("<INS:")
        ins_seq = record.INFO["SEQ"]
        record.ALT = [ref_seqs[record.CHROM][record.POS] + ins_seq]
        record.REF = record.ALT[0][0]


def _fix_aggregated_sv_calls(ref_fasta, infile, outfile):
    """Takes VCF made by graphtyper sv, which has big indels reported over
    multiple lines using the INS and DEL tags in the ALT column. The first
    line should be "AGGREGATED", then the next ones are other evidence. We want
    to just keep the AGGREGATED one, turning it into a 'normal' VCF line with
    the REF and ALT alleles having the indel, instead of in SEQ=... key/value in
    INFO fields"""
    found_aggregated = set()
    ref_seqs = {}
    for seq in pyfastaq.sequences.file_reader(ref_fasta):
        seq.id = seq.id.split()[0]
        ref_seqs[seq.id] = seq.seq

    with open(infile) as f_in, open(outfile, "w") as f_out:
        for line in f_in:
            if line.startswith("#"):
                print(line, end="", file=f_out)
                continue

            record = vcf_record.VcfRecord(line)

            # The SVs come in lines of 3, or as single records. When the SV is
            # 3 lines, we want to use the 'AGGREGATED' one. Otherwise we obviously
            # use the only record there is. When there's an aggregated, looks
            # like it is the first of the three - always check this.
            if record.ALT[0].startswith("<DEL:") or record.ALT[0].startswith("<INS:"):
                if "AGGREGATED" in record.ALT[0]:
                    assert record.ID not in found_aggregated
                    found_aggregated.add(record.ID)
                    _fix_sv_call(record, ref_seqs)
                else:
                    aggregated_id = record.ID.rsplit(".", maxsplit=1)[0]
                    if aggregated_id not in found_aggregated:
                        _fix_sv_call(record, ref_seqs)
                        found_aggregated.add(aggregated_id)
                    else:
                        continue

            print(record, file=f_out)


def run(outdir, ref_fasta, bam_file, vcfs, call_sv=False):
    outdir = os.path.abspath(outdir)
    bam_file = os.path.abspath(bam_file)
    original_dir = os.getcwd()
    vcfs = [os.path.abspath(x) for x in vcfs]
    os.mkdir(outdir)
    os.chdir(outdir)
    time_cmd = "/usr/bin/time -v"
    err_files = {}

    # Need to make a single VCF file as input to graphtyper. See
    # https://github.com/DecodeGenetics/graphtyper/issues/39
    vcf_for_graphtyper = "00.concat.vcf.gz"
    command = f"{time_cmd} bcftools concat -a -O v {' '.join(vcfs)} | bcftools sort -Oz -o {vcf_for_graphtyper}"
    utils.syscall(command, stdouterr="00.concat")
    err_files["bcftools_concat"] = "00.concat.err"
    utils.syscall(
        f"{time_cmd} tabix -p vcf {vcf_for_graphtyper}",
        stdouterr=f"{vcf_for_graphtyper}.tabix",
    )
    err_files["bcftools_concat_tabix"] = f"{vcf_for_graphtyper}.tabix.err"

    with open(ref_fasta) as f:
        ref_names = [x.rstrip().split()[0][1:] for x in f if x.startswith(">")]

    region_file = "00.regions.txt"
    with open(region_file, "w") as f:
        print(*ref_names, sep="\n", file=f)

    genotype_dir = "01.genotype"
    if call_sv:
        command = f"{time_cmd} graphtyper genotype_sv --output {genotype_dir} -vverbose --sam {bam_file} --region_file {region_file} --threads 1 {ref_fasta} {vcf_for_graphtyper}"
    else:
        command = f"{time_cmd} graphtyper genotype --output {genotype_dir} -vverbose --sam {bam_file} --region_file {region_file} --threads 1 {ref_fasta} --vcf {vcf_for_graphtyper}"
    utils.syscall(command, stdouterr=genotype_dir)
    err_files["genotype"] = f"{genotype_dir}.err"

    # graphtyper puts the results of each "region", which in our case is a
    # contig, into multiple VCFs file. Need to combine them into one.
    # graphtyper docs say to use "bcftools concat --naive".
    # Howwever, bcftools --naive option forces output to be compressed.
    # So don't use that # option!
    vcfs_to_concat = []
    for ref_name in ref_names:
        vcfs = glob.glob(os.path.join(genotype_dir, ref_name, "*.vcf.gz"))
        vcfs_to_concat.extend(vcfs)

    merged_vcf = "02.merged.vcf"
    command = f"{time_cmd} bcftools concat -a -O v {' '.join(vcfs_to_concat)} | bcftools sort -o {merged_vcf}"
    utils.syscall(command, stdouterr="02.bcftools_merge")
    err_files["bcftools_merge_2"] = "02.bcftools_merge.err"
    utils.time_and_memory_from_multiple_files(
        err_files, "resources.breakdown.json", "resources.json"
    )

    final_vcf = "02.final.vcf"
    _fix_aggregated_sv_calls(ref_fasta, merged_vcf, final_vcf)
    os.chdir(original_dir)


def run_combi(outdir, ref_fasta, bam_file, vcfs):
    outdir = os.path.abspath(outdir)
    bam_file = os.path.abspath(bam_file)
    original_dir = os.getcwd()
    vcfs = [os.path.abspath(x) for x in vcfs]
    os.mkdir(outdir)
    os.chdir(outdir)
    time_cmd = "/usr/bin/time -v"
    err_files = {}
    vcfs_svs = []
    vcfs_non_svs = []
    for vcf in vcfs:
        vcfs_svs.append(f"vcfs.{len(vcfs_svs)}.sv.vcf")
        vcfs_non_svs.append(f"vcfs.{len(vcfs_non_svs)}.non_sv.vcf")
        utils.split_vcf_into_svs_and_not_svs(
            vcf, vcfs_svs[-1], vcfs_non_svs[-1], len_cutoff=50
        )
        utils.bgzip_and_tabix_vcf(vcfs_svs[-1])
        utils.bgzip_and_tabix_vcf(vcfs_non_svs[-1])
        vcfs_svs[-1] += ".gz"
        vcfs_non_svs[-1] += ".gz"

    sv_dir = "call_svs"
    non_sv_dir = "call_non_svs"
    run(sv_dir, ref_fasta, bam_file, vcfs_svs, call_sv=True)
    run(non_sv_dir, ref_fasta, bam_file, vcfs_non_svs, call_sv=False)
    sv_vcf = os.path.join(sv_dir, "02.final.vcf")
    non_sv_vcf = os.path.join(non_sv_dir, "02.final.vcf")
    final_vcf = "final.vcf"
    utils.combine_sv_and_non_sv_vcfs(sv_vcf, non_sv_vcf, final_vcf)
    with open(os.path.join(sv_dir, "resources.json")) as f:
        resources = json.load(f)

    with open(os.path.join(non_sv_dir, "resources.json")) as f:
        new_resources = json.load(f)
    assert sorted(list(resources.keys())) == sorted(list(new_resources.keys()))
    for k, v in new_resources.items():
        resources[k] += v
    with open("resources.json", "w") as f:
        json.dump(resources, f, indent=2, sort_keys=True)
    os.chdir(original_dir)

