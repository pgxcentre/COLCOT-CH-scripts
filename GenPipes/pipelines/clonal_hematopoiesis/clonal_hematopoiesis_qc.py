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
from pipelines.clonal_hematopoiesis import clonal_hematopoiesis


log = logging.getLogger(__name__)

class ClonalHematopoiesisQC(clonal_hematopoiesis.ClonalHematopoiesis):
    """
    Clonal Hematopoiesis Amplicon Pipeline QC
    =================
    The Clonal Hematopoiesis Pipeline is based on the DNA-Seq pipeline and identifies somatic variants
    within target sequencing data. Here we compile QC metrics of the general pipeline to detect problems
    related to indexes, lab manipulation and more.


    """

    def __init__(self, protocol=None):
        self._protocol=protocol
        # Add pipeline specific arguments
        super(ClonalHematopoiesisQC, self).__init__(protocol)


    def umi_extract_qc(self):
        """
         1. After fgbio ExtractUmisFromBam to Extract read counts per design and per base in design intervals
        """
        jobs = []
        design_intervals_bed = config.param('amplicon_metrics', 'design_intervals', required=False)
		# run metrics  from readset names
        for readset in self.readsets:
            metrics_dir = os.path.join("metrics", "dna", readset.sample.name)
            picard_directory = os.path.join(metrics_dir, "picard_metrics")
            alignment_directory = os.path.join("alignment", readset.sample.name, readset.name)
            alignment_file_prefix = os.path.join(alignment_directory, readset.name)
            unmapped_bam = alignment_file_prefix + ".unmapped.bam"
            mapped_umi_unfiltered_bam =   alignment_file_prefix + ".merged.sorted.unfiltered.umi.bam"
            unmapped_umi_bam = alignment_file_prefix + ".unmapped.umi.bam"
            readset_bam = alignment_file_prefix + ".merged.sorted.umi.bam"
            output = alignment_file_prefix + ".mapped.umi.bam"
            output_samtools_amplicons = os.path.join(picard_directory, re.sub("\.bam$", ".depth.tsv", os.path.basename(output)))
            if design_intervals_bed == "" :
                design_intervals_bed = bvatools.resolve_readset_coverage_bed(readset)
            # Create jobs
            #jobs.append(
            #    concat_jobs([
            #        bash.mkdir(picard_directory),
            #        bvatools.depth_of_coverage(
            #            output,
            #            os.path.join(picard_directory, readset.sample.name + ".raw_coverage.tsv"),
            #            bvatools.resolve_readset_coverage_bed(readset),
            #            other_options=config.param('bvatools_depth_of_coverage', 'other_options', required=False)
            #            )
            #        ],
            #        name="bvatools_depth_of_coverage." + readset.name,
            #        samples=[readset.sample.name]
            #        )
            #    )
            jobs.append(
                concat_jobs([
                    bash.mkdir(picard_directory),
                    samtools.depth(
                        output,
                        output_samtools_amplicons,
                        bvatools.resolve_readset_coverage_bed(readset),
			"-J -aa"
                        )
                    ],
                    name="samtools_depth." + readset.name,
                    samples=[readset.sample.name]
                    )
                )

            #jobs.append(
            #    concat_jobs([
            #        bash.mkdir(picard_directory),
            #        bedtools.coverage(
            #            output,
            #            output_bedtools,
            #            design_intervals_bed,
            #            )
            #        ],
            #        name="amplicon_metrics." + readset.name,
            #        samples=[readset.sample.name]
            #        )
            #    )

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
            input_fastqc = [ os.path.join(metrics_directory, sample.name, "fastqc", readset) for readset in samples[sample.name]]
            input_dep += input_fastqc
            inputs += [os.path.join(metrics_directory, sample.name)]

        output = os.path.join(metrics_directory, "multiqc_report")

        job = multiqc.run(
            inputs,
            output,
            input_dep
            )
        job.name = "multiqc_all_samples"
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
            self.umi_extract_qc,
            self.run_multiqc
        ]

if __name__ == '__main__':
    argv = sys.argv
    if '--wrap' in argv:
        utils.utils.container_wrapper_argparse(argv)
    else:
        ClonalHematopoiesisQC()
