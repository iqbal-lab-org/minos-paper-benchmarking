import json
import logging
import operator
import os
import subprocess
import sys

from cluster_vcf_records import vcf_record, vcf_file_read
import pyfastaq


def touch(filename):
    with open(filename, "w") as f:
        pass


def rm_rf(filename):
    subprocess.check_output(f"rm -rf {filename}", shell=True)


def syscall(command, stdouterr=None):
    logging.info(f"Run command: {command}")
    completed_process = subprocess.run(
        command,
        shell=True,
        stderr=subprocess.PIPE,
        stdout=subprocess.PIPE,
        universal_newlines=True,
    )
    logging.info(f"Return code: {completed_process.returncode}")
    if completed_process.returncode != 0:
        print("Error running this command:", command, file=sys.stderr)
        print("cwd:", os.getcwd(), file=sys.stderr)
        print("Return code:", completed_process.returncode, file=sys.stderr)
        print(
            "Output from stdout:", completed_process.stdout, sep="\n", file=sys.stderr
        )
        print(
            "Output from stderr:", completed_process.stderr, sep="\n", file=sys.stderr
        )
        raise RuntimeError("Error in system call. Cannot continue")

    logging.info(f"stdout:\n{completed_process.stdout.rstrip()}")
    logging.info(f"stderr:\n{completed_process.stderr.rstrip()}")

    if stdouterr is not None:
        with open(f"{stdouterr}.out", "w") as f:
            print(completed_process.stdout.rstrip(), file=f)
        with open(f"{stdouterr}.err", "w") as f:
            print(completed_process.stderr.rstrip(), file=f)

    return completed_process


def bash_out_to_time_and_memory(infile):
    """Parses output file `out` of command run like this:
          /usr/bin/time -v foo &> out
    returns a dictionary of time and memory used"""
    stats = {}

    with open(infile) as f:
        in_time_lines = False
        for line in f:
            if not in_time_lines:
                if line.startswith("\tCommand being timed:"):
                    in_time_lines = True
            elif line.startswith("\tElapsed (wall clock) time (h:mm:ss or m:ss): "):
                time_fields = line.rstrip().split()[-1].split(":")
                if not 2 <= len(time_fields) <= 3:
                    raise Exception(f"Error getting time from this line: {line}")
                time_in_seconds = float(time_fields[-1]) + 60 * float(time_fields[-2])
                if len(time_fields) == 3:
                    time_in_seconds += 60 * 60 * float(time_fields[0])
                stats["wall_clock_time"] = time_in_seconds
            elif line.startswith("\tUser time (seconds): "):
                stats["user_time"] = float(line.rstrip().split()[-1])
            elif line.startswith("\tSystem time (seconds): "):
                stats["system_time"] = float(line.rstrip().split()[-1])
            elif line.startswith("\tMaximum resident set size (kbytes): "):
                stats["ram"] = float(line.rstrip().split()[-1])

    return stats


def json_to_file(dictionary, filename):
    with open(filename, "w") as f:
        json.dump(dictionary, f, indent=2, sort_keys=2)


def load_json(filename):
    with open(filename) as f:
        return json.load(f)


def time_and_memory_from_multiple_files(filenames, out_per_file, out_summary):
    summary = {}
    per_file = {k: bash_out_to_time_and_memory(v) for k, v in filenames.items()}
    for key in ("wall_clock_time", "user_time", "system_time"):
        summary[key] = sum([per_file[stage][key] for stage in per_file])
    summary["ram"] = max([per_file[stage]["ram"] for stage in per_file])
    json_to_file(per_file, out_per_file)
    json_to_file(summary, out_summary)


def time_and_memory_to_json(infile, outfile):
    stats = bash_out_to_time_and_memory(infile)
    json_to_file(stats, outfile)


def bgzip_and_tabix_vcf(infile, outfile=None):
    if outfile is None:
        outfile = infile + ".gz"
    syscall(
        f"/usr/bin/time -v bgzip -c {infile} > {outfile}", stdouterr=f"{outfile}.bgzip"
    )
    syscall(f"/usr/bin/time -v tabix -p vcf {outfile}", stdouterr=f"{outfile}.tabix")


def mean_depth_and_length_from_minos_vcf(infile):
    """Gets mean depth and max read length from header of VCF file made by minos"""
    max_length = None
    mean_depth = None

    with open(infile) as f:
        for line in f:
            if line.startswith("##minosMeanReadDepth"):
                mean_depth = float(line.rstrip().split("=")[-1])
            elif line.startswith("##minos_max_read_length"):
                max_length = float(line.rstrip().split("=")[-1])
            elif line.startswith("#CHROM"):
                break

    if max_length is not None and mean_depth is not None:
        return mean_depth, max_length
    else:
        raise RuntimeError(
            f"Error getting minos mean read depth and/or length from VCF file {infile}"
        )


# The next two functions are for for bayestyper.
# It requires the contig lines in the header of all
# input VCF files are the same, and then the records in the VCF are sorted
# in that order. bcftools and Picard can't do this, so have to do our own
# implementation. (Note you might think Picard can, but I couldn't get
# it to work because SortVcf assumes that there is at least one variant call
# for each contig, which is not always true.)
def fasta_to_ordered_names_and_lengths(fasta_in):
    refs_in_order = []
    ref_lengths = {}
    seq_reader = pyfastaq.sequences.file_reader(fasta_in)
    for seq in seq_reader:
        name = seq.id.split()[0]
        assert name not in ref_lengths
        refs_in_order.append(name)
        ref_lengths[name] = len(seq)
    return refs_in_order, ref_lengths


def sort_vcf_file(vcf_in, vcf_out, ref_lengths, ref_names):
    vcf_header = []
    vcf_records = {}
    with open(vcf_in) as f:
        for line in f:
            if line.startswith("#"):
                if line.startswith("##contig="):
                    continue
                vcf_header.append(line)
            else:
                chrom, pos, the_rest = line.split("\t", maxsplit=2)
                pos = int(pos)
                if chrom not in vcf_records:
                    vcf_records[chrom] = []
                vcf_records[chrom].append((pos, the_rest))
    assert vcf_header[-1].startswith("#CHROM")

    with open(vcf_out, "w") as f:
        for line in vcf_header[:-1]:
            print(line, end="", file=f)
        for name in ref_names:
            print(f"##contig=<ID={name},length={ref_lengths[name]}>", file=f)
        print(vcf_header[-1], end="", file=f)

        for name in ref_names:
            if name in vcf_records:
                vcf_records[name].sort()
                for (pos, the_rest) in vcf_records[name]:
                    print(name, pos, the_rest, sep="\t", end="", file=f)


def split_vcf_into_svs_and_not_svs(infile, out_sv, out_non_sv, len_cutoff=50):
    with open(infile) as f_in, open(out_sv, "w") as f_out_sv, open(
        out_non_sv, "w"
    ) as f_out_non_sv:
        for line in f_in:
            if line.startswith("#"):
                print(line, end="", file=f_out_sv)
                print(line, end="", file=f_out_non_sv)
            else:
                record = vcf_record.VcfRecord(line)
                length_diffs = [abs(len(record.REF) - len(x)) for x in record.ALT]
                if max(length_diffs) >= len_cutoff:
                    print(line, end="", file=f_out_sv)
                else:
                    print(line, end="", file=f_out_non_sv)


def combine_sv_and_non_sv_vcfs(sv_vcf, non_sv_vcf, outfile):
    pos_in_svs = {}
    variants = {}
    header = []
    with open(sv_vcf) as f:
        for line in f:
            if line.startswith("#"):
                header.append(line.rstrip())
            else:
                record = vcf_record.VcfRecord(line)
                # `graphtyper genotype_sv` can output lots of null calls that
                # do get called by `graphtyper genotype`. So remove the null
                # ones to get the best results. This should be rare on real
                # data, but def happens on the test data here. And does no
                # harm to remove the ./. calls.
                if record.FORMAT.get("GT", "") == "./.":
                    continue
                if record.CHROM not in variants:
                    variants[record.CHROM] = []
                    pos_in_svs[record.CHROM] = set()
                variants[record.CHROM].append(record)
                for i in range(record.POS, record.ref_end_pos() + 1, 1):
                    pos_in_svs[record.CHROM].add(i)

    with open(non_sv_vcf) as f:
        for line in f:
            if not line.startswith("#"):
                record = vcf_record.VcfRecord(line)
                use_record = True
                if record.CHROM in pos_in_svs:
                    for i in range(record.POS, record.ref_end_pos() + 1, 1):
                        if i in pos_in_svs[record.CHROM]:
                            use_record = False
                            break

                if use_record:
                    if record.CHROM not in variants:
                        variants[record.CHROM] = []
                    variants[record.CHROM].append(record)

    with open(outfile, "w") as f:
        print(*header, sep="\n", file=f)
        for chrom in variants:
            variants[chrom].sort(key=operator.attrgetter("POS"))
            for record in variants[chrom]:
                print(record, file=f)
