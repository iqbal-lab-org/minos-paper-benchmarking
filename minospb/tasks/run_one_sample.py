from minospb import sample_dir


def run(options):
    sdir = sample_dir.SampleDir(options.sample_dir)
    sdir.run_var_callers(
        options.reads1,
        options.reads2,
        options.ref_dir,
        options.sample_name,
        cortex_mem_height=22,
        cpus=1,
    )
    sdir.run_mergers(
        options.ref_dir,
        options.sample_name,
        kmc_ram=options.ram,
        minos_only=options.clockwork,
        minos_ref_splits=options.minos_ref_splits,
    )

    if options.clockwork:
        return

    sdir.eval_calls(
        options.ref_dir,
        options.truth_fasta,
        truth_vcf=options.truth_vcf,
        ref_mask_bed=options.ref_mask_bed,
        truth_mask_bed=options.truth_mask_bed,
        varifier_no_maxmatch=options.varifier_no_maxmatch,
    )
