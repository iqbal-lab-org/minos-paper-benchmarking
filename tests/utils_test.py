import filecmp
import os
import pytest

from minospb import utils

this_dir = os.path.dirname(os.path.abspath(__file__))
data_dir = os.path.join(this_dir, "data", "utils")


def test_fasta_to_ordered_names_and_lengths():
    infile = os.path.join(data_dir, "fasta_to_ordered_names_and_lengths.fa")
    expect_names = ["seq1", "seq2"]
    expect_lengths = {"seq1": 1, "seq2": 2}
    got_names, got_lengths = utils.fasta_to_ordered_names_and_lengths(infile)
    assert got_names == expect_names
    assert got_lengths == expect_lengths


def test_sort_vcf_file():
    infile = os.path.join(data_dir, "sort_vcf_file.in.vcf")
    expect_vcf = os.path.join(data_dir, "sort_vcf_file.expect.vcf")
    tmp_out = "tmp.sort_vcf_file.vcf"
    utils.rm_rf(tmp_out)
    ref_lengths = {"ref1": 42, "ref2": 100, "ref3": 11}
    ref_names = ["ref3", "ref2", "ref1"]
    utils.sort_vcf_file(infile, tmp_out, ref_lengths, ref_names)
    assert filecmp.cmp(tmp_out, expect_vcf, shallow=False)
    os.unlink(tmp_out)


def test_split_vcf_into_svs_and_not_svs():
    infile = os.path.join(data_dir, "split_vcf_into_svs_and_not_svs.in.vcf")
    tmp_out_sv = "tmp.split_vcf_into_svs_and_not_svs.sv.vcf"
    tmp_out_non_sv = "tmp.split_vcf_into_svs_and_not_svs.non_sv.vcf"
    utils.rm_rf(tmp_out_sv)
    utils.rm_rf(tmp_out_non_sv)
    utils.split_vcf_into_svs_and_not_svs(
        infile, tmp_out_sv, tmp_out_non_sv, len_cutoff=3
    )
    expect_sv = os.path.join(data_dir, "split_vcf_into_svs_and_not_svs.out_sv.vcf")
    expect_non_sv = os.path.join(
        data_dir, "split_vcf_into_svs_and_not_svs.out_non_sv.vcf"
    )
    assert filecmp.cmp(tmp_out_sv, expect_sv, shallow=False)
    assert filecmp.cmp(tmp_out_non_sv, expect_non_sv, shallow=False)
    os.unlink(tmp_out_sv)
    os.unlink(tmp_out_non_sv)


def test_combine_sv_and_non_sv_vcfs():
    sv_vcf = os.path.join(data_dir, "combine_sv_and_non_sv_vcfs.in.sv.vcf")
    non_sv_vcf = os.path.join(data_dir, "combine_sv_and_non_sv_vcfs.in.non_sv.vcf")
    expect_vcf = os.path.join(data_dir, "combine_sv_and_non_sv_vcfs.expect.vcf")
    tmp_out = "tmp.combine_sv_and_non_sv_vcfs.vcf"
    utils.rm_rf(tmp_out)
    utils.combine_sv_and_non_sv_vcfs(sv_vcf, non_sv_vcf, tmp_out)
    assert filecmp.cmp(tmp_out, expect_vcf, shallow=False)
    os.unlink(tmp_out)
