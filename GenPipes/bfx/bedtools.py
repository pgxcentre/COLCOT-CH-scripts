#!/usr/bin/env python

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
import os

# MUGQIC Modules
from core.config import *
from core.job import *

from bfx import ucsc

def graph(input_bam, output_bed_graph, library_type="PAIRED_END"):

    if library_type == "PAIRED_END":
        if "forward" in output_bed_graph:
            samtools_options = "-F 256 -f 81 "
        elif "reverse" in output_bed_graph:
            samtools_options = "-F 256 -f 97 "
        else:
            raise Exception("Error: PAIRED_END library was provided but no strand orientation could be determined from " + output_bed_graph + "...")
    else:
        samtools_options = "-F 256"

    return Job(
        [input_bam],
        [output_bed_graph],
        [
            ['bedtools', 'module_samtools'],
            ['bedtools', 'module_bedtools']
        ],
        command="""\
nmblines=$(samtools view {samtools_options} {input_bam} | wc -l) && \\
scalefactor=0$(echo "scale=2; 1 / ($nmblines / 10000000);" | bc) && \\
genomeCoverageBed {other_options} -bg -split -scale $scalefactor \\
  -ibam {input_bam} \\
  -g {chromosome_size} \\
  > {output_bed_graph}""".format(
            samtools_options=samtools_options,
            input_bam=input_bam,
            chromosome_size=config.param('bedtools_graph', 'chromosome_size', type='filepath'),
            other_options=config.param('bedtools_graph', 'other_options', required=False),
            output_bed_graph=output_bed_graph
        )
    )

def intersect(input_bam, output_bam, target_bed, include_header=False):

    return Job(
        [input_bam],
        [output_bam],
        [
            ['bedtools', 'module_bedtools']
        ],
        command="""\
bedtools intersect {other_options} {include_header} \\
  -a {input_bam} \\
  -b {target_bed} \\
  > {output_bam}""".format(
            input_bam=input_bam,
            target_bed=target_bed,
            other_options=config.param('bedtools_intersect', 'other_options', required=False),
            include_header="-header" if include_header else "",
            output_bam=output_bam
        )
    )

def intersect_beds(bed1, bed2, output_bed, other_options=""):

    return Job(
        [bed1, bed2],
        [output_bed],
        [
            ['bedtools', 'module_bedtools']
        ],
        command="""\
bedtools intersect {other_options} \\
  -a {bed1} \\
  -b {bed2}{output_bed}""".format(
            bed1 = bed1,
            bed2 = bed2,
            other_options = other_options,
            output_bed=" \\\n  > " + output_bed if output_bed else ""
        )
    )

def bamtobed(input_bam, output_bed):

    return Job(
        [input_bam],
        [output_bed],
        [
            ['bedtools', 'module_bedtools']
        ],
        command="""\
bedtools bamtobed {other_options} \\
  -i {input_bam}{output_bed}""".format(
            input_bam=input_bam,
            other_options=config.param('bedtools_coverage', 'other_options', required=False),
            output_bed=" \\\n  > " + output_bed if output_bed else ""
        )
    )

def bamtobedpe(input_bam, output_bed, other_options):

    return Job(
        [input_bam],
        [output_bed],
        [
            ['bedtools', 'module_bedtools']
        ],
        command="""\
bedtools bamtobed {other_options} \\
  -i {input_bam}{output_bed}""".format(
            input_bam=input_bam,
            other_options=other_options,
            output_bed=" \\\n  > " + output_bed if output_bed else ""
        )
    )

def coverage(input_file, output_file, interval=None):
    if interval is None :
        interval=config.param('bedtools_coverage', 'gc_intervals')

    return Job(
        [input_file, interval],
        [output_file],
        [
            ['bedtools', 'module_bedtools']
        ],
        command="""\
bedtools coverage {other_options} \\
  -a {intervals} \\
  -b {input} \\
  > {output_file}""".format(
            intervals=interval,
            input=input_file,
            other_options=config.param('bedtools_coverage', 'other_options', required=False),
            output_file=output_file
        )
    )

def genomecov(input_file, output_file):

    return Job(
        [input_file],
        [output_file],
        [
            ['bedtools', 'module_bedtools']
        ],
        command="""\
bedtools genomecov {other_options} \\
  {input} \\
  {genome} \\
  > {output_file}""".format(
            input="-ibam " + input_file if re.search("\.bam$", os.path.basename(input_file)) else "-i " + input_file,
            genome="-g " + config.param('bedtools_genomecov', 'genome_fasta') if not re.search("\.bam$", os.path.basename(input_file)) else "",
            other_options=config.param('bedtools_genomecov', 'other_options', required=False),
            output_file=output_file
        )
    )

def merge(input_bed, output_bed, other_options):

    return Job(
        [input_bed],
        [output_bed],
        [
            ['bedtools', 'module_bedtools']
        ],
        command="""\
bedtools merge \\
  -i {input_bed} {other_options}{output_bed}""".format(
            input_bed=input_bed,
            other_options=other_options,
            output_bed=" \\\n  > " + output_bed if output_bed else ""
        )
    )

