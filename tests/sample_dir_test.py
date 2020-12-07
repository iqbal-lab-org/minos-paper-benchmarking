import os
import pytest

from minospb import sample_dir, utils

this_dir = os.path.dirname(os.path.abspath(__file__))
data_dir = os.path.join(this_dir, "data", "sample_dir")


def test_sample_dir():
    outdir = "tmp.test_sample_dir"
    utils.rm_rf(outdir)
    reads1 = os.path.join(data_dir, "reads_1.fq")
    reads2 = os.path.join(data_dir, "reads_2.fq")
    ref_dir = os.path.join(data_dir, "Ref")
    truth_fasta = os.path.join(data_dir, "truth_ref.fa")
    ref_mask = os.path.join(data_dir, "ref_mask.bed")
    sample_name = "sample_42"
    sdir = sample_dir.SampleDir(outdir)

    assert not sdir.is_done("var_call")
    sdir.run_var_callers(
        reads1, reads2, ref_dir, sample_name, cortex_mem_height=20, cpus=1
    )
    assert sdir.is_done("var_call")
    sdir.run_var_callers(
        reads1, reads2, ref_dir, sample_name, cortex_mem_height=20, cpus=1
    )
    assert sdir.is_done("var_call")

    assert not sdir.is_done("merge")
    sdir.run_mergers(ref_dir, sample_name, kmc_ram=1)
    assert sdir.is_done("merge")
    sdir.run_mergers(ref_dir, sample_name, kmc_ram=1)
    assert sdir.is_done("merge")
    for filename in sdir.final_vcfs.values():
        assert os.path.exists(os.path.join(outdir, filename))

    assert not sdir.is_done("eval_calls")
    sdir.eval_calls(ref_dir, truth_fasta, ref_mask_bed=ref_mask)
    assert sdir.is_done("eval_calls")
    sdir.eval_calls(ref_dir, truth_fasta)
    assert sdir.is_done("eval_calls")

    assert os.path.exists(sdir.resources_json)
    assert os.path.exists(sdir.results_summary_json)

    utils.rm_rf(outdir)
