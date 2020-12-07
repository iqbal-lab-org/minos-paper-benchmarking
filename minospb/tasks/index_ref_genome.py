from clockwork import reference_dir


def run(options):
    ref_dir = reference_dir.ReferenceDir(directory=options.outdir)
    ref_dir.make_index_files(
        options.fasta_in,
        False,  # genome is not big
        True,  # make cortex index
        cortex_mem_height=22,
    )
