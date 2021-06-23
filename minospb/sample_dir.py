import glob
import logging
import os

from clockwork import cortex, read_map, read_trim
from cluster_vcf_records import vcf_record
import pyfastaq

from minospb import bayestyper, graphtyper, utils


class SampleDir:
    def __init__(self, root_dir):
        self.root_dir = os.path.abspath(root_dir)
        self.var_call_dir = "Varcall"
        self.merge_dir = "Merge_calls"
        self.eval_calls_dir = "Eval_calls"
        logging.info(f"Initialising sample directory: {self.root_dir}")

        if not os.path.exists(self.root_dir):
            logging.info(f"Make directory: {self.root_dir}")
            os.mkdir(self.root_dir)

        for d in self.var_call_dir, self.merge_dir, self.eval_calls_dir:
            full_path = os.path.join(self.root_dir, d)
            if not os.path.exists(full_path):
                logging.info(f"Make directory: {full_path}")
                os.mkdir(full_path)

        self.prefixes = {
            "trim_reads": "trim_reads",
            "map_reads": "map_reads",
            "samtools": os.path.join(self.var_call_dir, "samtools"),
            "cortex": os.path.join(self.var_call_dir, "cortex"),
            "bayestyper": os.path.join(self.merge_dir, "bayestyper"),
            "graphtyper_default": os.path.join(self.merge_dir, "graphtyper_default"),
            "graphtyper_sv": os.path.join(self.merge_dir, "graphtyper_sv"),
            "minos": os.path.join(self.merge_dir, "minos"),
        }

        self.completed_file = os.path.join(self.root_dir, "stages_complete")
        if os.path.exists(self.completed_file):
            with open(self.completed_file) as f:
                self.stages_complete = set([x.rstrip() for x in f])
        else:
            self.stages_complete = set()

        self.time_cmd = "/usr/bin/time -v"
        self.trim_reads1 = "trim_reads.1.fq.gz"
        self.trim_reads2 = "trim_reads.2.fq.gz"
        self.rmdup_bam = self.prefixes["map_reads"] + ".rmdup.bam"
        self.samtools_vcf = self.prefixes["samtools"] + ".vcf"
        self.samtools_vcf_gz = self.samtools_vcf + ".gz"
        self.cortex_dir = self.prefixes["cortex"]
        self.cortex_vcf = self.prefixes["cortex"] + ".vcf"
        self.cortex_vcf_gz = self.cortex_vcf + ".gz"
        self.minos_dir = self.prefixes["minos"]
        self.resources_json = os.path.join(self.root_dir, "resources.json")
        self.results_summary_json = os.path.join(self.root_dir, "results_summary.json")

        self.final_vcfs = {
            "bayestyper": os.path.join(self.prefixes["bayestyper"], "05.final.vcf"),
            "cortex": self.cortex_vcf,
            "graphtyper_default": os.path.join(
                self.prefixes["graphtyper_default"], "02.final.vcf"
            ),
            "graphtyper_sv": os.path.join(
                self.prefixes["graphtyper_sv"], "02.final.vcf"
            ),
            "minos": os.path.join(self.minos_dir, "final.vcf"),
            "samtools": self.samtools_vcf,
        }
        logging.info(
            f"Finish initialising sample directory. Stages done: {sorted(list(self.stages_complete))}"
        )

    def set_done(self, stage):
        self.stages_complete.add(stage)
        with open(self.completed_file, "w") as f:
            print(*sorted(list(self.stages_complete)), sep="\n", file=f)
        logging.info(f"Finished: {stage}")

    def is_done(self, stage):
        done = stage in self.stages_complete
        logging.info(f"is_done: {stage} {done}")
        return done

    def fix_cortex_vcf(self, original_cortex_vcf, ref_fasta):
        """The cortex VCF needs the "MISMAPPED_UNPLACEABLE" records removed
        because they can break other tools.
        It also needs any records removing where POS is bigger than the
        length of the ref, or is <1"""
        ref_seqs = {}
        pyfastaq.tasks.file_to_dict(ref_fasta, ref_seqs)
        ref_seqs = {k.split()[0]: v for k, v in ref_seqs.items()}

        with open(original_cortex_vcf) as f_in, open(self.cortex_vcf, "w") as f_out:
            for line in f_in:
                if line.startswith("#CHROM"):
                    for name, seq in sorted(ref_seqs.items()):
                        print(f"##contig=<ID={name},length={len(seq)}>", file=f_out)
                elif not line.startswith("#"):
                    try:
                        record = vcf_record.VcfRecord(line)
                    except:
                        continue

                    if (
                        "MISMAPPED_UNPLACEABLE" in record.FILTER
                        or record.POS < 1
                        or record.POS > len(ref_seqs[record.CHROM]) - 1
                        or not record.ref_string_matches_ref_sequence(
                            ref_seqs[record.CHROM]
                        )
                    ):
                        continue

                print(line, end="", file=f_out)

    def run_var_callers(
        self, reads1, reads2, ref_dir, sample_name, cortex_mem_height=22, cpus=1,
    ):
        if self.is_done("var_call"):
            return

        reads1 = os.path.abspath(reads1)
        reads2 = os.path.abspath(reads2)
        ref_dir = os.path.abspath(ref_dir)
        ref_fasta = os.path.join(ref_dir, "ref.fa")
        original_dir = os.getcwd()
        os.chdir(self.root_dir)

        if not self.is_done("trim_reads"):
            logging.info("Start trim_reads")
            read_trim.run_trimmomatic(
                reads1,
                reads2,
                self.trim_reads1,
                self.trim_reads2,
                qual_trim="LEADING:10 TRAILING:10 SLIDINGWINDOW:4:15",
            )
            self.set_done("trim_reads")

        if not self.is_done("map_reads"):
            logging.info("Start map_reads")
            read_map.map_reads(
                ref_fasta,
                self.trim_reads1,
                self.trim_reads2,
                self.rmdup_bam,
                rmdup=True,
                read_group=("1", sample_name),
            )
            utils.syscall(f"samtools index {self.rmdup_bam}")
            self.set_done("map_reads")

        ref_names, ref_lengths = utils.fasta_to_ordered_names_and_lengths(ref_fasta)

        # --------------------------- samtools -------------------------------------
        if not self.is_done("var_call.samtools"):
            logging.info("Start var_call.samtools")
            utils.rm_rf(self.prefixes["samtools"] + "*")
            command = f"samtools mpileup -ugf {ref_fasta} {self.rmdup_bam} | bcftools call -vm -O v -o {self.samtools_vcf}"
            utils.syscall(command, stdouterr=self.prefixes["samtools"])
            utils.sort_vcf_file(
                self.samtools_vcf, self.samtools_vcf, ref_lengths, ref_names
            )
            utils.bgzip_and_tabix_vcf(self.samtools_vcf)
            self.set_done("var_call.samtools")

        # ---------------------------- cortex --------------------------------------
        if not self.is_done("var_call.cortex"):
            logging.info("Start var_call.cortex")
            utils.rm_rf(self.prefixes["cortex"] + "*")
            ctx = cortex.CortexRunCalls(
                ref_dir,
                self.rmdup_bam,
                self.cortex_dir,
                sample_name,
                mem_height=cortex_mem_height,
            )
            ctx.run()
            vcfs = [
                x
                for x in glob.glob(
                    os.path.join(self.cortex_dir, "cortex.out", "vcfs", "**raw.vcf")
                )
            ]
            assert len(vcfs) == 1
            self.fix_cortex_vcf(vcfs[0], ref_fasta)
            utils.sort_vcf_file(
                self.cortex_vcf, self.cortex_vcf, ref_lengths, ref_names
            )
            utils.bgzip_and_tabix_vcf(self.cortex_vcf)
            self.set_done("var_call.cortex")

        self.set_done("var_call")
        os.chdir(original_dir)

    def run_mergers(
        self, ref_dir, sample_name, minos_only=False, minos_ref_splits=1, kmc_ram=12
    ):
        logging.info("Start run_mergers")
        assert self.is_done("var_call")
        if self.is_done("merge"):
            return

        resource_dirs = {}
        ref_dir = os.path.abspath(ref_dir)
        ref_fasta = os.path.join(ref_dir, "ref.fa")
        original_dir = os.getcwd()
        os.chdir(self.root_dir)

        # --------------------------- bayestyper -----------------------------------
        if not minos_only and not self.is_done("merge.bayestyper"):
            logging.info("Start merge.bayestyper")
            utils.rm_rf(self.prefixes["bayestyper"] + "*")
            bayestyper.run(
                self.prefixes["bayestyper"],
                ref_fasta,
                self.rmdup_bam,
                self.samtools_vcf,
                self.cortex_vcf,
                sample_name,
                kmc_ram=kmc_ram,
            )
            resource_dirs["bayestyper"] = self.prefixes["bayestyper"]
            self.set_done("merge.bayestyper")

        # ----------------------------- minos --------------------------------------
        if not self.is_done("merge.minos"):
            logging.info("Start merge.minos")
            utils.rm_rf(self.prefixes["minos"] + "*")
            split_opt = (
                "" if minos_ref_splits == 1 else f"--total_splits {minos_ref_splits}"
            )
            command = f"{self.time_cmd} minos --debug adjudicate {split_opt} --force --reads {self.rmdup_bam} {self.minos_dir} {ref_fasta} {self.samtools_vcf} {self.cortex_vcf}"
            utils.syscall(command, stdouterr=self.prefixes["minos"])
            json_out = os.path.join(self.minos_dir, "resources.json")
            resource_dirs["minos"] = self.minos_dir
            utils.time_and_memory_to_json(f"{self.prefixes['minos']}.err", json_out)
            old = self.prefixes["minos"]
            new = os.path.join(self.minos_dir, "stdouterr")
            os.rename(f"{old}.out", f"{new}.out")
            os.rename(f"{old}.err", f"{new}.err")
            self.set_done("merge.minos")

        # --------------------------- graphtyper -----------------------------------
        if not minos_only and not self.is_done("merge.graphtyper"):
            logging.info("Start merge.graphtyper_default")
            utils.rm_rf(self.prefixes["graphtyper_default"] + "*")
            graphtyper.run(
                self.prefixes["graphtyper_default"],
                ref_fasta,
                self.rmdup_bam,
                [self.samtools_vcf_gz, self.cortex_vcf_gz],
            )
            resource_dirs["graphtyper_default"] = self.prefixes["graphtyper_default"]
            self.set_done("merge.graphtyper_default")

        # --------------------------- graphtyper_sv --------------------------------
        if not minos_only and not self.is_done("merge.graphtyper_sv"):
            logging.info("Start merge.graphtyper_sv")
            utils.rm_rf(self.prefixes["graphtyper_sv"] + "*")
            graphtyper.run(
                self.prefixes["graphtyper_sv"],
                ref_fasta,
                self.rmdup_bam,
                [self.samtools_vcf_gz, self.cortex_vcf_gz],
                call_sv=True,
            )
            resource_dirs["graphtyper_sv"] = self.prefixes["graphtyper_sv"]
            self.set_done("merge.graphtyper_sv")

        resources = {
            k: utils.load_json(os.path.join(v, "resources.json"))
            for k, v in resource_dirs.items()
        }
        utils.json_to_file(resources, self.resources_json)

        self.set_done("merge")
        os.chdir(original_dir)

    def eval_calls(
        self,
        ref_dir,
        truth_fasta,
        truth_vcf=None,
        ref_mask_bed=None,
        truth_mask_bed=None,
        varifier_no_maxmatch=False,
    ):
        logging.info("Start eval_calls")
        assert self.is_done("var_call")
        assert self.is_done("merge")
        if self.is_done("eval_calls"):
            return

        varifier_opts = []
        if truth_vcf is not None:
            varifier_opts.extend(["--truth_vcf", os.path.abspath(truth_vcf)])
        if varifier_no_maxmatch:
            varifier_opts.append("--no_maxmatch")

        varifier_opts = " ".join(varifier_opts)

        ref_fasta = os.path.abspath(os.path.join(ref_dir, "ref.fa"))
        truth_fasta = os.path.abspath(truth_fasta)
        results = {}
        mask_opt_string = ""
        if ref_mask_bed is not None:
            ref_mask_bed = os.path.abspath(ref_mask_bed)
            mask_opt_string = f"--ref_mask {ref_mask_bed}"
        if truth_mask_bed is not None:
            truth_mask_bed = os.path.abspath(truth_mask_bed)
            mask_opt_string += f" --truth_mask {truth_mask_bed}"
        original_dir = os.getcwd()
        os.chdir(self.root_dir)

        for tool_name, vcf_to_eval in self.final_vcfs.items():
            tooldir = os.path.join(self.eval_calls_dir, tool_name)
            if not os.path.exists(tooldir):
                os.mkdir(tooldir)
            results[tool_name] = {}

            for filter_opt in "no_filter", "PASS":
                results[tool_name][filter_opt] = {}
                for mask in ["masked", "unmasked"]:
                    outdir = os.path.join(tooldir, f"{mask}.{filter_opt}")
                    done = f"eval_calls.{tool_name}.{mask}.{filter_opt}"
                    logging.info(f"Start {done}")
                    summary_json = os.path.join(outdir, "summary_stats.json")
                    if not self.is_done(done):
                        if mask == "masked" and mask_opt_string == "":
                            utils.rm_rf(outdir)
                            os.mkdir(outdir)
                            utils.json_to_file({}, summary_json)
                        else:
                            filter_pass = (
                                ""
                                if filter_opt == "no_filter"
                                else f"--filter_pass {filter_opt}"
                            )
                            mask_opt = "" if mask == "unmasked" else mask_opt_string
                            command = f"varifier vcf_eval {varifier_opts} {mask_opt} {filter_pass} --force {truth_fasta} {ref_fasta} {vcf_to_eval} {outdir}"
                            utils.syscall(command)
                        self.set_done(done)

                    results[tool_name][filter_opt][mask] = utils.load_json(summary_json)

        utils.json_to_file(results, self.results_summary_json)
        self.set_done("eval_calls")
        os.chdir(original_dir)
