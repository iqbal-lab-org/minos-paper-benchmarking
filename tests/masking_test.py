import filecmp
import os
import pytest


from minospb import masking, utils

this_dir = os.path.dirname(os.path.abspath(__file__))
data_dir = os.path.join(this_dir, "data", "masking")


def test_pysam_pileup_list_is_ok():
    assert not masking.pysam_pileup_list_is_ok([], "A", 1, 100.0)
    assert masking.pysam_pileup_list_is_ok(["A"], "A", 1, 100.0)
    assert masking.pysam_pileup_list_is_ok(["A"], "a", 1, 100.0)
    assert not masking.pysam_pileup_list_is_ok(["A"], "n", 1, 100.0)
    assert not masking.pysam_pileup_list_is_ok(["A"], "N", 1, 100.0)
    assert not masking.pysam_pileup_list_is_ok(["A"], "A", 2, 100.0)
    assert not masking.pysam_pileup_list_is_ok(["A"], "C", 1, 100.0)
    assert masking.pysam_pileup_list_is_ok(["a"], "A", 1, 100.0)
    assert not masking.pysam_pileup_list_is_ok(["a"], "C", 1, 100.0)
    pileup_list = ["A"] * 4 + ["a"] * 4 + ["*", "a+4ACGT"]
    assert masking.pysam_pileup_list_is_ok(pileup_list, "A", 8, 80.0)
    assert not masking.pysam_pileup_list_is_ok(pileup_list, "A", 8, 80.1)
    assert not masking.pysam_pileup_list_is_ok(pileup_list, "A", 9, 80.0)

def test_bam_to_mask():
    ref_fa = os.path.join(data_dir, "ref.fa")
    bam_file = os.path.join(data_dir, "mapped_reads.bam")
    outfile = "tmp.masking.bam_to_mask.bed"
    utils.rm_rf(outfile)
    masking.bam_to_mask(bam_file, ref_fa, outfile, 5, 50)
    expect = os.path.join(data_dir, "bam_to_mask.percent50.bed")
    assert filecmp.cmp(outfile, expect, shallow=False)
    os.unlink(outfile)

    masking.bam_to_mask(bam_file, ref_fa, outfile, 5, 90)
    expect = os.path.join(data_dir, "bam_to_mask.percent90.bed")
    assert filecmp.cmp(outfile, expect, shallow=False)
    os.unlink(outfile)


def test_trim_and_map_reads_and_make_mask():
    ref_fa = os.path.join(data_dir, "ref.fa")
    reads1 = os.path.join(data_dir, "reads_1.fq.gz")
    reads2 = os.path.join(data_dir, "reads_2.fq.gz")
    outdir = "tmp.trim_map_mask"
    utils.rm_rf(outdir)
    got = masking.trim_and_map_reads_and_make_mask(
        ref_fa, reads1, reads2, outdir, 5, 90
    )
    expect = os.path.join(data_dir, "trim_and_map_reads_and_make_mask.bed")
    assert filecmp.cmp(got, expect, shallow=False)
    utils.rm_rf(outdir)
