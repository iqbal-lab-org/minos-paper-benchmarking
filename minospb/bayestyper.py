import os

from minospb import utils


# We're running on bacteria and Plasmodium. By default
# "BayesTyper assumes human ploidy levels.". This writes the input ploidy file
# to tell bayesTyper that all contigs are haploid.
def write_ploidy_file(ref_fasta, outfile):
    with open(ref_fasta) as f_in, open(outfile, "w") as f_out:
        for line in f_in:
            if line.startswith(">"):
                name = line.lstrip(">").rstrip().split()[0]
                print(name, 1, 1, sep="\t", file=f_out)


def run(outdir, ref_fasta, bam_file, vcfs, sample_name, kmc_ram=12):
    outdir = os.path.abspath(outdir)
    ref_fasta = os.path.abspath(ref_fasta)
    bam_file = os.path.abspath(bam_file)
    vcfs = [os.path.abspath(x) for x in vcfs]
    vcfs_norm_out = [f"02.typertools.{i}.norm" for i in range(len(vcfs))]
    vcfs_norm_vcfs = [f"{x}.vcf" for x in vcfs_norm_out]
    original_dir = os.getcwd()
    os.mkdir(outdir)
    os.chdir(outdir)
    time_cmd = "/usr/bin/time -v"
    kmc_out = "01.kmc"
    typertools_combine_out = "02.typertools.combine"
    typertools_cluster_out = "03.typertools.cluster"
    typertools_ploidy_file = "04.typertools.ploidy.tsv"
    typertools_genotype_out = "04.typertools.genotype"
    tsv_file = "03.typertools.samples.tsv"
    err_files = {}

    command = f"{time_cmd} kmc -t1 -m{kmc_ram} -k55 -ci1 -fbam {bam_file} {kmc_out} ."
    utils.syscall(command, stdouterr=kmc_out)
    err_files["kmc"] = f"{kmc_out}.err"

    command = f"{time_cmd} bayesTyperTools makeBloom -k {kmc_out}"
    utils.syscall(command, stdouterr="makeBloom")
    err_files["makeBloom"] = "makeBloom.err"
    vcf_opt_string = []

    for i, vcf_in in enumerate(vcfs):
        command = f"{time_cmd} bcftools norm -f {ref_fasta} {vcf_in} -o {vcfs_norm_vcfs[i]}"
        utils.syscall(command, stdouterr=vcfs_norm_out[i])
        err_files[f"norm_{i}"] = f"{vcfs_norm_out[i]}.err"
        vcf_opt_string.append(f"{i}:{vcfs_norm_vcfs[i]}")
    vcf_opt_string = ",".join(vcf_opt_string)
    command = f"{time_cmd} bayesTyperTools combine -v {vcf_opt_string} -o {typertools_combine_out} -z"
    utils.syscall(command, stdouterr=typertools_combine_out)
    err_files["combine"] = f"{typertools_combine_out}.err"

    with open(tsv_file, "w") as f:
        print(sample_name, "M", os.path.join(outdir, kmc_out), sep="\t", file=f)

    command = f"{time_cmd} bayesTyper cluster -r 42 -v {typertools_combine_out}.vcf.gz -s {tsv_file} -g {ref_fasta}"
    utils.syscall(command, stdouterr=typertools_cluster_out)
    err_files["cluster"] = f"{typertools_cluster_out}.err"

    write_ploidy_file(ref_fasta, typertools_ploidy_file)

    command = f"{time_cmd} bayesTyper genotype --noise-genotyping -r 42 -v bayestyper_unit_1/variant_clusters.bin -y {typertools_ploidy_file} -c bayestyper_cluster_data -s {tsv_file} -g {ref_fasta} -o bayestyper_unit_1/bayestyper"
    utils.syscall(command, stdouterr=typertools_genotype_out)
    err_files["genotype"] = f"{typertools_genotype_out}.err"

    utils.time_and_memory_from_multiple_files(
        err_files, "resources.breakdown.json", "resources.json"
    )

    final_vcf = os.path.join("bayestyper_unit_1", "bayestyper.vcf")
    assert os.path.exists(final_vcf)
    os.symlink(final_vcf, "05.final.vcf")
    os.chdir(original_dir)
