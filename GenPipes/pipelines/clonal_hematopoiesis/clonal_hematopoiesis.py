#!/usr/bin/env python
# -*- coding: utf-8 -*-

################################################################################
# Copyright (C) 2014, 2015 GenAP, McGill University and Genome Quebec Innovation Centre
#
# This file is part of MUGQIC Pipelines.
#
# MUGQIC Pipelines is free software: you can redistribute it and/or modify
# it under the terms of the GNU Lesser General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# MUGQIC Pipelines is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public License
# along with MUGQIC Pipelines.  If not, see <http://www.gnu.org/licenses/>.
################################################################################

# Python Standard Modules
import logging
import math
import os
import re
import sys
import collections

# Append mugqic_pipelines directory to Python library path
sys.path.append(os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(sys.argv[0])))))

# MUGQIC Modules
from core.config import config
from core.job import Job, concat_jobs, pipe_jobs

from bfx import adapters
from bfx import bcftools
from bfx import bwa
from bfx import bvatools
from bfx import bedtools
from bfx import cutadapt
from bfx import fastqc
from bfx import gatk
from bfx import gatk4
from bfx import igvtools
from bfx import multiqc
from bfx import picard
from bfx import picard2
from bfx import fgbio
from bfx import samtools
from bfx import sambamba
from bfx import tools
from bfx import vardict
from bfx import htslib
from bfx import vt
from bfx import snpeff
from bfx import vep
from bfx import bash_cmd as bash

from pipelines.dnaseq import dnaseq

log = logging.getLogger(__name__)

class ClonalHematopoiesis(dnaseq.DnaSeq):
    """
    Clonal Hematopoiesis Amplicon Pipeline
    =================
    The Clonal Hematopoiesis Pipeline is based on the DNA-Seq pipeline and identifies somatic variants
    within target sequencing data. Somatic variants are identified by allowing for varying ploidy in
    the form of allele fractions for each variant. If Unique Molecular Identifiers are detected the
    fgbio CallMolecularConsensusReads tool is used. Recalibration is not run either because typically,
    these datasets are targeted with amplicons or custom capture which render recalibration useless.

    Data Analysis Recommendations for Illumina Paired-end Sequencing using CleanPlex UMI Panels follows:

    1- Trim trailing adapters in the R1 and R2 reads using cutadapt
        cutadapt -g CCTACACGACGCTCTTCCGATCT \\
        -a AGATCGGAAGAGCACACGTCTGAA \\
        -A AGATCGGAAGAGCGTCGTGTAGG \\
        -G TTCAGACGTGTGCTCTTCCGATCT -o SampleID.1.fastq -p SampleID.2.fastq SampleID.rawreads.1.fastq SampleID.rawreads.2.fastq -e 0.1 -O 9 -m 20 -n 2
        Only when the panel is GCP53, cutadapt needs to be run for a second time:
        cutadapt -g GTTATTGGACGTTTA \\
        -a TTACCTTTTATGAAG \\
        -A TAAACGTCCAATAAC \\
        -G CTTCATAAAAGGTAA \\
        -e 0.1 -O 9 -m 20 -n 2 \\
        -o SampleID.trimmed.1.fastq -p SampleID.trimmed.1.fastq SampleID.1.fastq  SampleID.2.fastq

    2- Step 1: Convert the trimmed fastq reads to unmapped BAM
        java -jar picard.jar FastqToSam F1=SampleID.trimmed.1.fastq F2=SampleID.trimmed.2.fastq \
        SM=SampleID O=SampleID.step1-unmapped.bam

    3- Step 2: Extract the UMI information to store as a tag
        java -jar fgbio.jar ExtractUmisFromBam --input=SampleID.step1-unmapped.bam --output=SampleID.step2-unmapped-umi.bam --read-structure=16M+T 16M+T --single-tag=RX --molecular-index-tags=ZA ZB

    4- Step 3: Convert the template read to fastq
        java -jar picard.jar SamToFastq I=SampleID.step2-unmapped-umi.bam F=SampleID.step3-R1Template.fastq F2=SampleID.step3-R2Template.fastq

    5- Step 4: Map the template reads against the human genome
        bwa mem -t No_threads_computation  -M hg19.fa SampleID.step3-R1Template.fastq SampleID.step3-R2Template.fastq > SampleID.step4-MappedTemplate.bam

    6- Step 5: Merge the mapped reads with the UMI info
        java -jar picard.jar MergeBamAlignment UNMAPPED=SampleID.step2-unmapped-umi.bam ALIGNED=SampleID.step4-MappedTemplate.bam R=hg19.fa O=SampleID.step5-Mapped-UMI.bam

    7- Step 6: Group reads by UMI
        java -jar fgbio.jar GroupReadsByUmi --input=SampleID.step5-Mapped-UMI.bam --output=SampleID.step6-Mapped-GroupedByUMI.bam --strategy=adjacency --edits=1 --min-map=10

    8- Step 7: Call consensus on grouped reads
        java -jar fgbio.jar CallMolecularConsensusReads --input=SampleID.step6-Mapped-GroupedByUMI.bam --output=SampleID.step7-ConsensusReads-Unmapped.bam --min-reads=3

    9- Step 8: Convert the consensus reads to fastq
        java -jar picard.jar SamToFastq I=SampleID.step7-ConsensusReads-Unmapped.bam F=SampleID.step8-R1Consensus.fastq F2=SampleID.step8-R2Consensus.fastq

    10- Step 9: Map the consensus reads against the human genome
        bwa mem -t No_threads_computation -M hg19.fa SampleID.step8-R1Consensus.fastq SampleID.step8-R2Consensus.fastq > SampleID.step9-Mapped-Consensus.bam

    11- Step 10: Merge the unmapped and the mapped consensus files
        java -jar picard.jar MergeBamAlignment UNMAPPED=SampleID.step7-ConsensusReads-Unmapped.bam ALIGNED=SampleID.step9-Mapped-Consensus.bam R=$Ref O=SampleID.step10-Mapped-Consensus-UMI.bam

    12- Step 11: Create index
        samtools index -b SampleID.step10-Mapped-Consensus-UMI.bam SampleID.step10-Mapped-Consensus-UMI.bai

    13- Step 12: Call variant caller
        VarDict -G hg19.fa -N SampleID -b SampleID.step10-Mapped-Consensus-UMI.bam -f Allele_Freq -q 11 -r 2 -k 1 -I 50 -c 1 -S 2 -E 3 -g 4 PanelName.ampInsert.bed | teststrandbias.R | var2vcf_valid.pl -N $tumor_group -E -f $af > SampleID_AF_consensus.vardict.vcf

    """

    def __init__(self, protocol=None):
        self._protocol=protocol
        # Add pipeline specific arguments
        super(ClonalHematopoiesis, self).__init__(protocol)

    def fastqc(self):
        """
        [FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/) provide quality 
        control checks on raw sequence data coming from high throughput sequencing pipelines.
        It provides a modular set of analyses including: per base sequence quality, per sequence 
        quality scores, per base sequence content, per base GC content, per sequence GC content, 
        per base N contenf, Sequence length distribitions, sequence duplication, onverrepresented
        sequences and kmer content metrics.
        """
        jobs = []
        for readset in self.readsets:
            raw_fastq_directory = os.path.join("raw_reads", readset.sample.name)
            raw_fastq_file_prefix = os.path.join(raw_fastq_directory, readset.name)
	        # Create jobs
            fastqc_directory = os.path.join("metrics", "dna", readset.sample.name, "fastqc")
            output_dir = os.path.join(fastqc_directory)
            if readset.run_type == "PAIRED_END":
                candidate_input_files = []
                if readset.fastq1 and readset.fastq2:
                    candidate_input_files.append([readset.fastq1, readset.fastq2])
                [fastq1, fastq2] = self.select_input_files(candidate_input_files)
                output = os.path.join(fastqc_directory, re.sub(".fastq.gz$", "_fastqc.zip", os.path.basename(fastq1)))
            elif readset.run_type == "SINGLE_END":
                candidate_input_files = []
                if readset.fastq1:
                    candidate_input_files.append([readset.fastq1])
                [fastq1] = self.select_input_files(candidate_input_files)
                fastq2 = None
                output = os.path.join(fastqc_directory, re.sub(".fastq.gz$", "_fastqc.zip", os.path.basename(fastq1)))

            else:
                _raise(SanitycheckError("Error: run type \"" + readset.run_type +
                "\" is invalid for readset \"" + readset.name + "\" (should be PAIRED_END or SINGLE_END)!"))
            # Create jobs
            adapter_file = config.param('fastqc', 'adapter_file', required=False, type='filepath')
            adapter_job = None

            if not adapter_file:
                adapter_file = os.path.join(output_dir, "adapter.tsv")
                adapter_job = adapters.create(
                    readset,
                    adapter_file,
                    fastqc=True
                )
            jobs.append(
                concat_jobs([
                    bash.mkdir(
                        output_dir,
                        remove=True
                        ),
                    adapter_job,
                    fastqc.fastqc(
                        fastq1,
                        fastq2,
                        output_dir,
                        output,
                        adapter_file
                        )
                    ],
                    name="fastqc." + readset.name,
                    samples=[readset.sample.name]
                    )
            )
        return jobs

    def cutadapt(self):
        """
        Raw reads quality trimming and removing adapters is performed using [Cutadapt](https://cutadapt.readthedocs.io/en/stable/index.html).
        'Adapter1' and 'Adapter2' columns from the readset file ar given to Cutadapt. For PAIRED_END readsets, both adapters are used.
        For SINGLE_END readsets, only Adapter1 is used and left unchanged.
        To trim the front of the read use adapter_5p_fwd and adapter_5p_rev (for PE only) in cutadapt section of ini file.

        This step takes as input files:

        1. FASTQ files from the readset file if available
        2. Else, FASTQ output files from previous picard_sam_to_fastq conversion of BAM files
        """
        jobs = []
        for readset in self.readsets:
            trim_directory = os.path.join("trim", readset.sample.name)
            trim_file_prefix = os.path.join(trim_directory, readset.name)
            metrics_directory = os.path.join("metrics", "dna", readset.sample.name, "cutadapt")
            if readset.run_type == "PAIRED_END":
                candidate_input_files = []
                if readset.fastq1 and readset.fastq2:
                    candidate_input_files.append([readset.fastq1, readset.fastq2])
                [fastq1, fastq2] = self.select_input_files(candidate_input_files)
                adapter_fwd = readset.adapter1
                adapter_rev = readset.adapter2
            elif readset.run_type == "SINGLE_END":
                candidate_input_files = []
                if readset.fastq1:
                    candidate_input_files.append([readset.fastq1])
                [fastq1] = self.select_input_files(candidate_input_files)
                fastq2 = None
                adapter_fwd = readset.adapter1
                adapter_rev = None
            else:
                _raise(SanitycheckError("Error: run type \"" + readset.run_type +
                "\" is invalid for readset \"" + readset.name + "\" (should be PAIRED_END or SINGLE_END)!"))
            jobs.append(
                concat_jobs([
                    bash.mkdir(trim_directory),
                    bash.mkdir(metrics_directory),
                    cutadapt.trim(
                        fastq1,
                        fastq2,
                        trim_file_prefix,
                        adapter_fwd,
                        adapter_rev,
                        metrics_directory
                        )
                    ],
                    name="cutadapt." + readset.name)
                )
        return jobs

    def metrics(self):
        """
        Compute metrics and generate coverage tracks per sample. Multiple metrics are computed at this stage:
        Number of raw reads, Number of filtered reads, Number of aligned reads, 
        Median, mean and standard deviation of insert sizes of reads after alignment, percentage of bases
        covered at X reads (%_bases_above_50 means the % of exons bases which have at least 50 reads)
        whole genome or targeted percentage of bases covered at X reads (%_bases_above_50 means the % of exons
        bases which have at least 50 reads). 
        """
        library = {}
        for readset in self.readsets:
            if not library.has_key(readset.sample):
                library[readset.sample] = "SINGLE_END"
            if readset.run_type == "PAIRED_END":
                library[readset.sample] = "PAIRED_END"

        jobs = []
        for sample in self.samples:
            input_file_prefix = os.path.join("alignment", sample.name, sample.name + ".sorted.")
            input = input_file_prefix + "bam"
            picard_directory = os.path.join("metrics", "dna", sample.name, "picard_metrics")
            mkdir_job = bash.mkdir(picard_directory, remove=True)
            output_bedtools = os.path.join(picard_directory, re.sub("\.bam$", ".BedGraph", os.path.basename(input)))
            # Create jobs
            jobs.append(
                concat_jobs([
                    mkdir_job,
                    picard2.collect_multiple_metrics(
                        input,
                        os.path.join(picard_directory, sample.name + ".all.metrics"),
                        library_type=library[sample]
                        ),
                    bvatools.depth_of_coverage(
                        input,
                        os.path.join(picard_directory, sample.name + ".coverage.tsv"),
                        bvatools.resolve_readset_coverage_bed(sample.readsets[0]),
                        other_options=config.param('bvatools_depth_of_coverage', 'other_options', required=False)
                        ),
                    picard2.collect_gcbias_metrics(input, 
                                                os.path.join(picard_directory, sample.name + ".qcbias_metrics.txt"),
                                                os.path.join(picard_directory, sample.name + ".qcbias_metrics.pdf"),
                                                os.path.join(picard_directory, sample.name + ".qcbias_summary_metrics.txt")),
                    bedtools.coverage(
                        input,
                        output_bedtools
                        )
                    ],
                    name="picard_collect_multiple_metrics." + sample.name,
                    samples=[sample]
                    )
                )

        return jobs

    def picard_calculate_hs_metrics(self):
        """
        Compute on target percent of hybridisation based capture.
        """
        jobs = []
        created_interval_lists = []

        for sample in self.samples:
            coverage_bed = bvatools.resolve_readset_coverage_bed(sample.readsets[0])
            if coverage_bed:
                metrics_directory = os.path.join("metrics", "dna", sample.name, "picard_metrics")
                interval_list = re.sub("\.[^.]+$", ".interval_list", os.path.basename(coverage_bed))
                if not interval_list in created_interval_lists:
                    job = tools.bed2interval_list(None, coverage_bed, interval_list)
                    job.name = "interval_list." + os.path.basename(coverage_bed)
                    jobs.append(job)
                    created_interval_lists.append(interval_list)
                # Create metrics directory to collect with multiqc
                mkdir_job = bash.mkdir(metrics_directory, remove=True)
                alignment_directory = os.path.join("alignment", sample.name)
                [input] = self.select_input_files([
                    [os.path.join(alignment_directory, sample.name + ".sorted.bam")]
                ])
                # Run hs metrics
                jobs.append(
                    concat_jobs([
                        mkdir_job,
                        picard2.collect_hs_metrics(
                            input,
                            os.path.join(metrics_directory, re.sub("bam$", "onTarget.tsv", os.path.basename(input))),
                            interval_list
                        )
                    ],
                    name="picard_calculate_hs_metrics." + sample.name,
                    samples=[sample]
                    )
                )
        return jobs

    def umi_extract(self):
        """
         1. Use fgbio ExtractUmisFromBam to Extract the UMI information and store as a tag
        """
        jobs = []
        for readset in self.readsets:
            trim_directory = os.path.join("trim", readset.sample.name)
            trim_file_prefix = os.path.join(trim_directory, readset.name)
            alignment_directory = os.path.join("alignment", readset.sample.name, readset.name)
            alignment_file_prefix = os.path.join(alignment_directory, readset.name)
            metrics_dir = os.path.join("metrics", "dna", readset.sample.name)
            picard_directory = os.path.join(metrics_dir, "picard_metrics")
            unmapped_bam = alignment_file_prefix + ".unmapped.bam"
            mapped_umi_unfiltered_bam =   alignment_file_prefix + ".merged.sorted.unfiltered.umi.bam"
            unmapped_umi_bam = alignment_file_prefix + ".unmapped.umi.bam"
            readset_bam = alignment_file_prefix + ".merged.sorted.umi.bam"
            output = alignment_file_prefix + ".mapped.umi.bam"
            if readset.run_type == "PAIRED_END":
                fastq = trim_file_prefix + ".trim.pair1.fastq.gz"
                second_end_fastq = trim_file_prefix + ".trim.pair2.fastq.gz"
                fastq_umi = trim_file_prefix + ".trim.umi.pair1.fastq.gz"
                second_end_fastq_umi = trim_file_prefix + ".trim.umi.pair2.fastq.gz"

            if readset.run_type == "SINGLE_END":
                fastq = trim_file_prefix + ".trim.single.fastq.gz"
                second_end_fastq = None
                fastq_umi = trim_file_prefix + ".trim.single.fastq.gz"
                second_end_fastq_umi = None
            #fastq_to_sam(fastq, second_end_fastq, output, sample_name=None, read_group=None, library=None, platform=None, run=None, center=None)
            # Create jobs
            jobs.append(
                concat_jobs([
                    bash.mkdir(alignment_directory),
                    picard2.fastq_to_sam(
                        fastq,
                        second_end_fastq,
                        unmapped_bam,
                        readset.sample.name,
                        readset.name,
                        readset.library if readset.library else readset.sample.name,
                        config.param('bwa_mem', 'sequencing_technology') if config.param('bwa_mem', 'sequencing_technology', required=False) else "Illumina",
                        readset.run + "_" + readset.lane if readset.run and readset.lane else "",
                        config.param('bwa_mem', 'sequencing_center') if config.param('bwa_mem', 'sequencing_center', required=False) else ""
                    ),
                    fgbio.extractumisfrombam(
                            unmapped_bam,
                            unmapped_umi_bam
                            ),
                    picard2.sam_to_fastq(
                        unmapped_umi_bam,
                        fastq_umi,
                        second_end_fastq_umi,
                        )],
                    name="picard_fastq_to_sam." + readset.name,
                    samples=[readset.sample.name]
                 )
            )
            jobs.append(
                concat_jobs([
                    pipe_jobs([
                        bwa.mem(
                            fastq_umi,
                            second_end_fastq_umi,
                            read_group="'@RG" + \
                                "\\tID:" + readset.name + \
                                "\\tSM:" + readset.sample.name + \
                                "\\tLB:" + (readset.library if readset.library else readset.sample.name) + \
                                ("\\tPU:run" + readset.run + "_" + readset.lane if readset.run and readset.lane else "") + \
                                ("\\tCN:" + config.param('bwa_mem', 'sequencing_center') if config.param('bwa_mem', 'sequencing_center', required=False) else "") + \
                                ("\\tPL:" + config.param('bwa_mem', 'sequencing_technology') if config.param('bwa_mem', 'sequencing_technology', required=False) else "Illumina") + \
                                "'"
                                ),
                        sambamba.view(
                            "/dev/stdin",
                            "/dev/stdout",
                            "-S -f bam"
                            ),
                        sambamba.sort(
                            "/dev/stdin",
                            readset_bam,
                            config.param('sambamba_sort_sam', 'tmp_dir', required=True)
                            )
                    ]),
                    picard2.merge_alignments(
                        readset_bam,
                        unmapped_umi_bam,
                        output)
		    #sambamba.view(
                    #        input_bam=mapped_umi_unfiltered_bam,
                    #        output_bam=output,
                    #        options=config.param('fgbio_groupreadsbyumi', 'options_filter_bam', required=True)
                    #        )
                    ],
                    name="bwa_mem." + readset.name,
                    samples=[readset.sample.name]
                    )
                )

        return jobs

    def umi_deduplicate(self):
        """
         1. Use fgbio GroupReadsByUmi to group reads that appeared to have come from the same original
         molecule (this groups reads by template as well as by UMI sequence.)
         2. Used fgbio CallDuplexConsensusReads to create duplex consensus reads by combining the 
         evidence from a top and bottom single stranded consensus read
        """

        jobs = []
        for readset in self.readsets:
            trim_directory = os.path.join("trim", readset.name)
            trim_file_prefix = os.path.join(trim_directory, readset.name)
            alignment_directory = os.path.join("alignment", readset.sample.name, readset.name)
            alignment_file_prefix = os.path.join(alignment_directory, readset.name)
            [input] = self.select_input_files([
                [alignment_file_prefix + ".mapped.umi.bam"],
            ])
            dedup_directory =  os.path.join("deduplication", readset.name)
            metrics_dir = os.path.join("metrics", "dna", readset.sample.name, "fgbio")
            fastq1=os.path.join(dedup_directory, readset.name + '.sort-umis_corrected-cumi-1.fastq.gz')
            fastq2=os.path.join(dedup_directory, readset.name + '.sort-umis_corrected-cumi-2.fastq.gz')
            bam_unmapped = alignment_file_prefix + ".consensus.unmapped.bam"
            bam_unmapped_sorted = alignment_file_prefix + ".consensus.unmapped_sorted.bam"
            bam_grouped = alignment_file_prefix + ".consensus.mapped.grouped.bam"
            bam_mapped = alignment_file_prefix + ".consensus.mapped.bam"
            output = alignment_file_prefix + ".sorted.bam"
            index_bam = alignment_file_prefix + ".sorted.bam.bai"
            jobs.append(
                concat_jobs([
                    bash.mkdir(dedup_directory),
                    bash.mkdir(metrics_dir),
                    fgbio.groupreadsbyumi(
                            input,
                            bam_grouped,
                            os.path.join(metrics_dir, readset.sample.name + '.stats.tsv')
                        ),
                    fgbio.callmolecularconsensusreads(
                           bam_grouped,
                           bam_unmapped,
                           readset.name
                        ),
					pipe_jobs([
                        samtools.view(
							bam_unmapped,
                            None,
                            "-bh -F 2048"
                            ),
                        fgbio.sortbam(
                            '/dev/stdin',
                            bam_unmapped_sorted
                            )
					]),
                    picard2.sam_to_fastq(
                        bam_unmapped_sorted,
                        fastq1,
                        fastq2
                        ),
                    pipe_jobs([
                        bwa.mem(
                            fastq1,
                            fastq2,
                            read_group="'@RG" + \
                                "\\tID:" + readset.name + \
                                "\\tSM:" + readset.sample.name + \
                                "\\tLB:" + (readset.library if readset.library else readset.sample.name) + \
                                #("\\tPU:" + readset.name) + \
                                ("\\tPU:" + readset.run + "." + readset.lane if readset.run and readset.lane else "") + \
                                ("\\tCN:" + config.param('bwa_mem', 'sequencing_center') if config.param('bwa_mem', 'sequencing_center', required=False) else "") + \
                                ("\\tPL:" + config.param('bwa_mem', 'sequencing_technology') if config.param('bwa_mem', 'sequencing_technology', required=False) else "Illumina") + \
                                "'"
                                ),
						samtools.view(
                            "-",
                            None,
                            "-bh -F 2048"
                            ),
						fgbio.sortbam(
                            None,
                            bam_mapped
                            )
                        ]),
                    picard2.merge_alignments(
                        bam_mapped,
                        bam_unmapped_sorted,
                        output),
                    sambamba.index(
                        output,
                        index_bam
                        )
                    ],
                    name="umi_deduplicate." + readset.name,
                    samples=[readset.sample.name]
                    )
                )
        # Return
        return jobs


    def filter_alignements(self):
        """
        Filter problematic reads using Sambamba. This can be caused by primer dimers
        Compute metrics and generate coverage tracks per sample. Multiple metrics are computed at this stage:
        Number of raw reads, Number of filtered reads, Number of aligned reads, 
        Median, mean and standard deviation of insert sizes of reads after alignment, percentage of bases
        covered at X reads (%_bases_above_50 means the % of exons bases which have at least 50 reads)
        whole genome or targeted percentage of bases covered at X reads (%_bases_above_50 means the % of exons
        bases which have at least 50 reads). 
        """
        jobs = []
        for sample in self.samples:
            input_file_prefix = os.path.join("alignment", sample.name, sample.name + ".sorted.")
            input = input_file_prefix + "bam"
            picard_directory = os.path.join("metrics", "dna", sample.name, "picard_metrics")
            mkdir_job = bash.mkdir(picard_directory, remove=True)
            output = os.path.join(picard_directory, re.sub("\.bam$", ".filtered.bam", os.path.basename(input)))
            index_bam = os.path.join(picard_directory, re.sub("\.bam$", ".filtered.bai", os.path.basename(input)))

            # Create jobs
            jobs.append(
                concat_jobs([
                    mkdir_job,
					sambamba.view(
                            input_bam=input,
                            output_bam=output,
                            options=config.param('post_filter_alignments', 'options_filter_bam', required=True)
                            ),
                    sambamba.index(
                        output,
                        index_bam
                        )
                    ],
                    name="post_filter_alignments." + sample.name,
                    samples=[sample.name]
                    )
                )
        return jobs

    def variants_mutect2(self):
        """
        GATK MuTect2 caller for SNVs and Indels, tumor only workflow:
        https://docs.gdc.cancer.gov/Data/Bioinformatics_Pipelines/DNA_Seq_Variant_Calling_Pipeline/#tumor-only-variant-calling-workflow.

        """
        jobs = []

        nb_jobs = config.param('gatk_mutect2', 'nb_jobs', type='posint')
        pon = config.param('gatk_mutect2', 'panel_of_normals', type='filepath')
        if nb_jobs > 50:
            log.warning("Number of mutect jobs is > 50. This is usually much. Anything beyond 20 can be problematic.")

        # Create interval lists
        created_interval_lists = []
        for sample in self.samples:
            tumor_alignment_directory = os.path.join("alignment", sample.name)
            mutect_directory = os.path.join("variants", sample.name, "rawMuTect2")
            input_tumor = self.select_input_files(
                [[os.path.join(tumor_alignment_directory, sample.name + ".sorted.dup.recal.bam")],
                 [os.path.join(tumor_alignment_directory, sample.name + ".sorted.dup.bam")],
                 [os.path.join(tumor_alignment_directory, sample.name + ".sorted.bam")]])
            interval_list = None
            # Sequencing_artifacts_metrics file
            output_prefix = os.path.join(mutect_directory, sample.name)
            # pileup summaries
            output_pileup = output_prefix + '.targeted_sequencing.table'
            # contamination 
            output_contamination = output_prefix + '.targeted_sequencing.contamination.table'
            # raw mutect output
            output_mutect2 = os.path.join(mutect_directory, sample.name + ".mutect2.vcf.gz")
            output_mutect2_stats = os.path.join(mutect_directory, sample.name + ".mutect2.vcf.gz.stats")
            # sorted vcf
            output_mutect2_sorted = os.path.join(mutect_directory, sample.name + ".mutect2.targeted_sequencing.mutect2.tumor_only.sorted.vcf.gz")
            # filtered vcf
            output_mutect2_filtered = os.path.join(mutect_directory, sample.name + ".mutect2.targeted_sequencing.mutect2.tumor_only.contFiltered.vcf.gz")
            # final vcf
            output_mutect2_final = os.path.join(mutect_directory, sample.name + ".mutect2.targeted_sequencing.mutect2.raw_somatic_mutation.vcf.gz")
            # Check if bed file is specified in readset table
            coverage_bed = bvatools.resolve_readset_coverage_bed(sample.readsets[0])

            if coverage_bed:
                interval_list = re.sub("\.[^.]+$", ".interval_list", coverage_bed)
                if not interval_list in created_interval_lists:
                    job = tools.bed2interval_list(
                        None,
                        coverage_bed,
                        interval_list
                    )
                    job.name = "interval_list." + os.path.basename(coverage_bed)
                    jobs.append(job)
                    created_interval_lists.append(interval_list)

            jobs.append(concat_jobs([
                # Create output directory since it is not done by default by GATK tools
                bash.mkdir(
                mutect_directory,
                remove=True
                ),
                gatk4.collect_sequencing_artifacts_metrics(
                    input_tumor[0],
                    output_prefix),
                gatk4.get_pileup_summaries(
                    input_tumor[0],
                    output_pileup,
                    interval_list),
                gatk4.calculate_contamination(
                    output_pileup,
                    output_contamination),
                gatk4.mutect2_tumor_only(
                    input_tumor[0],
                    sample.name,
                    output_mutect2,
                    panel_of_normals=pon,
                    interval_list=interval_list
                    ),
                gatk4.sort_vcfs(
                    [output_mutect2],
                    output_mutect2_sorted
                    )
              ], name="gatk_mutect2." + sample.name)
            )

            jobs.append(concat_jobs([
                gatk4.filter_mutect_calls(
                    output_mutect2_sorted,
                    output_mutect2_filtered,
                    contamination=output_contamination,
                    stats=output_mutect2_stats),
                gatk4.filter_by_orientation_bias(
                    output_mutect2_filtered,
                    output_mutect2_final,
                    output_prefix + '.bait_bias_summary_metrics')
                ], name="filter_mutect_calls." + sample.name)
            )
        return jobs

    def vardict(self):
        """
        Variant calls in single sample mode. [VarDict](https://github.com/AstraZeneca-NGS/VarDict) 
        implements amplicon bias aware variant calling from targeted sequencing experiments, rescue 
        of long indels by realigning bwa soft clipped reads and better scalability than many Java 
        based variant callers effect annotation. 
        """
        #somatic(input, sample_name, output=None, region=[]):
        #VarDict -G hg19.fa -N SampleID -b SampleID.step10-Mapped-Consensus-UMI.bam -f Allele_Freq -q 11 -r 2 -k 1 -I 50 -c 1 -S 2 -E 3 -g 4 PanelName.ampInsert.bed | 
        #teststrandbias.R |
        # var2vcf_valid.pl -N $tumor_group -E -f $af > SampleID_AF_consensus.vardict.vcf

        jobs = []
        bed_file = None
        bed_file = config.param('vardict', 'regions_bed', required=False)
        # command
        for sample in self.samples:
            tumor_alignment_directory = os.path.join("alignment", sample.name)
            vardict_directory = os.path.join("variants", sample.name, "vardict")
            input_tumor = self.select_input_files(
                [[os.path.join(tumor_alignment_directory, sample.name + ".sorted.dup.recal.bam")],
                 [os.path.join(tumor_alignment_directory, sample.name + ".sorted.dup.bam")],
                 [os.path.join(tumor_alignment_directory, sample.name + ".sorted.bam")]])
            interval_list = None
            if bed_file == "":
                bed_file = bvatools.resolve_readset_coverage_bed(sample.readsets[0])
            # Sequencing_artifacts_metrics file
            output = os.path.join(vardict_directory, sample.name + 'AF_consensus.vardict.vcf')
            output_index = output + '.gz'
            jobs.append(concat_jobs([
                        bash.mkdir(
                            vardict_directory,
                            remove=True
                        ),
                        pipe_jobs([
                            vardict.single(
                                input_tumor[0],
                                sample.name,
                                None,
                                bed_file
                            ),
                            Job([],
                                [],
                                [],
                                command="$VARDICT_BIN/teststrandbias.R",
                            ),
                            vardict.var2vcf_valid(
                                None,
                                sample.name,
                                output
                            )
                        ]),
                        htslib.bgzip_tabix(
                             output,
                             output_index
                            )

                    ], name="vardict." + sample.name ))

        return jobs

    def snp_effect(self):
        """
        Variant effect annotation. The .vcf files are annotated for variant effects using the SnpEff software.
        SnpEff annotates and predicts the effects of variants on genes (such as amino acid changes).
        """
        jobs = []
        output_directory = "variants"
        snpeff_prefix = os.path.join(output_directory, "allSamples")
        for sample in self.samples:
            vardict_directory = os.path.join("variants", sample.name, "vardict")
            # final vcf
            output_vardict_final = os.path.join(mutect_directory, sample.name + ".mutect2.targeted_sequencing.mutect2.raw_somatic_mutation.vcf.gz")
            snpeff_prefix = os.path.join(mutect_directory, sample.name + ".mutect2.targeted_sequencing.mutect2.raw_somatic_mutation_annotated")
            output_norm = snpeff_prefix + ".vt.vcf.gz"
            output_annotated  = snpeff_prefix + ".vt.snpeff.vcf"
            jobs.append(concat_jobs([
                Job(command="mkdir -p " + output_directory, samples=self.samples),
                pipe_jobs([
                    vt.decompose_and_normalize_mnps(
                        output_mutect2_final,
                        None
                    ),
                    htslib.bgzip_tabix(
                        None,
                        output_norm
                    )
                ]),
                snpeff.compute_effects(output_norm, output_annotated, split=True),
                htslib.bgzip_tabix( output_annotated, snpeff_prefix + ".gz")
                ], name="normalize_and_compute_effects." + sample.name)
            )
        return jobs


    def annotate_filter(self):
        """
        Variant effect annotation. The .vcf files are annotated for variant effects using the
        (Ensembl VEP software)[https://useast.ensembl.org/info/docs/tools/vep/script/index.html#contents].
        The Variant Effect Predictor (VEP) uses "cache files" to read genomic data for better performance.
        """
        jobs = []
        output_directory = "variants"
        if config.param('ensembl_vep', 'output_format', required=False):
        	output_format =  config.param('ensembl_vep', 'output_format', required=False)
        else:
        	output_format = 'vcf'
        options_split_vep = config.param('bcftools_split_vep', 'options', required=True)
        annotations_header = config.param('bcftools_split_vep', 'header', required=True)
        for sample in self.samples:
            vardict_directory = os.path.join("variants", sample.name, "vardict")
            metrics_dir = os.path.join("metrics", "dna", sample.name, "vep")
            # final vcf
            output_vardict_final = os.path.join(vardict_directory, sample.name + 'AF_consensus.vardict.vcf.gz')
            snpeff_prefix = os.path.join(metrics_dir, sample.name + 'AF_consensus.vardict')
            output_norm = snpeff_prefix + ".vt.vcf.gz"
            metrics_prefix = os.path.join(metrics_dir, sample.name)
            output_annotated = metrics_prefix + ".vep." + output_format + '.gz'
            output_annotated_tab = metrics_prefix + ".vep." + output_format + ".tsv"
            output_filtered = metrics_prefix + ".vep.filtered." + output_format

            jobs.append(concat_jobs([
                Job(command="mkdir -p " + metrics_dir, samples=self.samples),
                pipe_jobs([
                    vt.decompose_and_normalize_mnps(
                        output_vardict_final,
                        None
                    ),
                    htslib.bgzip_tabix(
                        None,
                        output_norm
                    )
                    ])
                ],
                name="decompose_and_normalize_mnps." + sample.name)
            )
            # annotate and filter
            jobs.append(concat_jobs([
                Job(command="mkdir -p " + metrics_dir, samples=self.samples),
                pipe_jobs([
                	vep.annotate(output_norm,
				'STDOUT'),
			        htslib.bgzip_tabix(
        	                None,
                	        output_annotated
	                    )
        	        ]),
                vep.filter(output_annotated, output_filtered),
		        bcftools.split_vep(output_annotated, output_annotated_tab, options_split_vep)
                ],
                name="ensembl_vep." + sample.name)
            )
        return jobs

    def amplicon_metrics(self):
        """
        Metrics used to check for the presence of non specific products in amplicons
        """
        jobs = []
        output_directory = "metrics"
        filter_bam_options =  config.param('amplicon_metrics', 'options_filter_bam', required=False)
        for sample in self.samples:
            metrics_directory = os.path.join("metrics", "dna", sample.name, "picard_metrics")
            alignment_directory = os.path.join("alignment", sample.name, sample.name)
            alignment_file_prefix = os.path.join(alignment_directory, sample.name)
            # final vcf
            output_alignment = alignment_file_prefix + ".mapped.umi.bam"
            output_bedpe = os.path.join(metrics_directory, sample.name + ".mapped.umi.f2only.sorted.fixed.bedpe")
            output_intersect = os.path.join(metrics_directory, sample.name + ".mapped.umi.f2only.sorted.fixed.bedpe.intersect.bed")
            output_intersect_inverse = os.path.join(metrics_directory, sample.name + ".mapped.umi.f2only.sorted.fixed.bedpe.intersect.outer.bed")

            jobs.append(concat_jobs([
                Job(command="mkdir -p " + metrics_directory, samples=self.samples),
                pipe_jobs([
		        sambamba.view(
                        output_alignment,
                        None,
                        '-f bam ' + filter_bam_options if filter_bam_options else  '-f bam '
                    ),
                    sambamba.sort(
                        "/dev/stdin",
                        "/dev/stdout",
                        config.param('amplicon_metrics', 'tmp_dir', required=True),
  			            config.param('amplicon_metrics', 'options_sort', required=True),
                        ),
                    samtools.fixmate(
                        "/dev/stdin",
                        None,
                        config.param('amplicon_metrics', 'options_fixmate')
                        ),
                    bedtools.bamtobedpe(
			            "stdin",
            			output_bedpe,
                        config.param('amplicon_metrics', 'options_bedtools', required=True)
		            )
                ]),
                ], name="amplicon_metrics." + sample.name)
            )
            jobs.append(concat_jobs([
                Job(command="mkdir -p " + metrics_directory, samples=self.samples),
                    pipe_jobs([
                        bash.awk(
                            output_bedpe,
                            None,
                            "'{if($1==$4) print $1 \"\\t\" $2 \"\\t\" $6 \"\\t\" $7}'"
                        ),
                        Job(command="sort -k1,1 -k2,2n "),
		                bedtools.intersect_beds(
                            "stdin",
                            bvatools.resolve_readset_coverage_bed(sample.readsets[0]),
	    	                None,
		                    other_options=config.param('amplicon_metrics', 'options_bed_intersect', required=True)
                        ),
                        Job(command="csvcut -t -H -c 5,6,7,8"),
                        Job(command="sort -k 4,4"),
                        Job(command="uniq -c"),
                        Job(output_files=[output_intersect],
                            module_entries=[ ['amplicon_metrics', 'module_python']],
                            command="csvformat -T" + " > " + output_intersect)
                    ]),
                ], name="amplicon_metrics.intersect." + sample.name)
            )
            jobs.append(concat_jobs([
                Job(command="mkdir -p " + metrics_directory, samples=self.samples),
                pipe_jobs([
                    bedtools.intersect_beds(
                        output_bedpe,
                        bvatools.resolve_readset_coverage_bed(sample.readsets[0]),
                        None,
                        other_options=config.param('amplicon_metrics', 'options_bed_intersect_inverse', required=True)
                    ),
                    Job(command="sort -k1,1 -k2,2n "),
                    bedtools.merge(
                        "stdin",
                        None,
                        "-c 1 -o count"
                    ),
                    Job(output_files=[output_intersect_inverse], command="sort -k4,4n > " + output_intersect_inverse)
                ]),
                ], name="amplicon_metrics.intersect.inverse." + sample.name)
            )
        return jobs

    def amplicon_metrics_novaseq(self):
        """
        Metrics used to get coverage per base metrics using samtools
        """
        jobs = []
        output_directory = "metrics"
        design_intervals_bed = config.param('amplicon_metrics', 'design_intervals', required=False)
        # samtools depth per base
        for sample in self.samples:
            metrics_directory = os.path.join("metrics", "dna", sample.name, "picard_metrics")
            alignment_directory = os.path.join("alignment", sample.name, sample.name)
            alignment_file_prefix = os.path.join(alignment_directory, sample.name)
            output_alignment = alignment_file_prefix + ".sorted.bam"
            if design_intervals_bed == "" :
                design_intervals_bed = bvatools.resolve_readset_coverage_bed(sample.readsets[0])
            # Create jobs
            output_samtools_amplicons = os.path.join(metrics_directory, re.sub("\.bam$", ".depth.tsv", os.path.basename(output_alignment)))
            jobs.append(concat_jobs([
                Job(command="mkdir -p " + metrics_directory, samples=self.samples),
                samtools.depth(
                    output_alignment,
                    output_samtools_amplicons,
                    design_intervals_bed,
                    "-J -aa"
                    )
                ], name="samtools_depth." + sample.name
			))

        return jobs

    def chip_filtering(self):
        """
        Generate whitelist and manual variants files
        """
        jobs = []
        output_directory = "metrics"
        module_python = config.param('chip_filtering', 'module_python', required=False)
        other_options = config.param('chip_filtering', 'other_options', required=False)
        reference = config.param('chip_filtering', 'genome_fasta', type='filepath')
        peptides = config.param('chip_filtering', 'peptides', type='filepath')

        # samtools depth per base
        for sample in self.samples:
            metrics_prefix = os.path.join("metrics", "dna", sample.name, "vep", sample.name)
            output_annotated_tab = metrics_prefix + ".vep.vcf.tsv"
            output_prefix =  output_annotated_tab + '.filtered'
            jobs.append(concat_jobs([
                Job(command="mkdir -p " + metrics_directory, samples=self.samples),
				Job([output_annotated_tab],
			        [output_prefix + 'whitelist.tsv'],
			        [module_python],
			        command="""\
clonal_hematopoiesis_variants_filter.py \\
-i {source} \\
-o {dest}\\
-g {reference}\\
-a {peptides}{other_options}""".format(
			        source=output_annotated_tab,
		            dest=output_prefix,
					reference=reference,
					peptides=peptides)
					)
                ], name="chip_filtering." + sample.name)
            )

        return jobs

    def run_multiqc(self):

        jobs = []

        metrics_directory = os.path.join("metrics", "dna")
        input_dep = []
        inputs = []
        samples = collections.defaultdict(list)
        for readset in self.readsets:
            samples[readset.sample.name].append(re.sub(".fastq.gz$", "_fastqc.zip", os.path.basename(readset.fastq1)))
            #samples[readset.sample.name].append(re.sub(".fastq.gz$", "_fastqc.zip", os.path.basename(readset.fastq2)))

        for sample in self.samples:
            input_qcbias = os.path.join(metrics_directory, sample.name, "picard_metrics", sample.name + ".qcbias_metrics.txt")
            input_cutadapt = os.path.join(metrics_directory, sample.name, "cutadapt", sample.name + ".log")
            input_all_picard = os.path.join(metrics_directory, sample.name, "picard_metrics", sample.name + ".all.metrics.alignment_summary_metrics")
            input_fastqc = [ os.path.join(metrics_directory, sample.name, "fastqc", readset) for readset in samples[sample.name]]

            input_dep += [
                input_qcbias,
                input_all_picard,
            ] + input_fastqc

            inputs += [os.path.join(metrics_directory, sample.name)]

        output = os.path.join(metrics_directory, "multiqc_report")

        job = multiqc.run(
            inputs,
            output,
            input_dep
            )
        job.name = "multiqc"
        job.samples = self.samples

        jobs.append(job)

        return jobs

    @property
    def steps(self):
        return [
            self.picard_sam_to_fastq,
            self.fastqc,
            self.cutadapt,
            self.umi_extract,
            self.umi_deduplicate,
            self.sambamba_merge_sam_files,
            self.cram_output,
            self.metrics,
            self.picard_calculate_hs_metrics,
            #self.variants_mutect2,
            #self.snp_effect,
            self.vardict,
            self.annotate_filter,
            self.run_multiqc,
  	        self.amplicon_metrics,
            self.amplicon_metrics_novaseq,
			self.filter_alignements
        ]

if __name__ == '__main__': 
    argv = sys.argv
    if '--wrap' in argv:
        utils.utils.container_wrapper_argparse(argv)
    else:
        ClonalHematopoiesis()
