from minospb import masking


def run(options):
    masking.trim_and_map_reads_and_make_mask(
        options.ref_fa,
        options.reads1,
        options.reads2,
        options.outdir,
        options.min_depth,
        options.min_depth_pc,
    )
