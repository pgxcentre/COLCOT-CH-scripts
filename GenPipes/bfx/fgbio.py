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
import logging
import os

# MUGQIC Modules
from core.config import *
from core.job import *

log = logging.getLogger(__name__)

def addumi(
    input_bam,
    input_umi,
    output_bam,
    output_bai
    ):

    inputs = [input_bam, input_umi]
    outputs = [output_bam,output_bai]
    return Job(
        inputs,
        outputs,
        [
            ['fgbio_addumi', 'module_java'],
            ['fgbio_addumi', 'module_fgbio']
        ],

        command="""\
java -Djava.io.tmpdir={tmp_dir} {java_other_options} -Xmx{ram} -jar $FGBIO_JAR AnnotateBamWithUmis \\
  --input {input_bam} \\
  --fastq {input_umi} \\
  --output {output_bam} \\
  {other_options}""".format(
        tmp_dir=config.param('fgbio_addumi', 'tmp_dir'),
        java_other_options=config.param('fgbio_addumi', 'java_other_options'),
        ram=config.param('fgbio_addumi', 'ram'),
        input_bam=input_bam,
        input_umi=input_umi,
        output_bam=output_bam,
	other_options=config.param('fgbio_addumi', 'other_options') if config.param('fgbio_addumi', 'other_options',required=False) else ""
        ),
        removable_files=[output_bam]
    )

def correct_readname(
    input_umi,
    output_umi_corrected
    ):

    inputs = [input_umi]
    outputs = [output_umi_corrected]

    if input_umi.lower().endswith('.gz') :
        input_opener="zcat "
    else:
        input_opener="cat "

    return Job(
        inputs,
        outputs,
        command="""\
{input_opener} {input_umi} | tr ' ' '_' | gzip -c - > {output_umi_corrected}""".format(
        input_opener=input_opener,
        input_umi=input_umi,
        output_umi_corrected=output_umi_corrected
        ),
        removable_files=[output_umi_corrected]
    )

def trim_primers(input_bam, output_bam, hard_clip=False):
    inputs = [input_bam]
    outputs = [output_bam]

    return Job(
        inputs,
        outputs,
        [
            ['fgbio_trim_primers', 'module_java'],
            ['fgbio_trim_primers', 'module_fgbio']
        ],
        command="""\
java -Djava.io.tmpdir={tmp_dir} {java_other_options} -Xmx{ram} -jar $FGBIO_JAR TrimPrimers \\
  --input {input_bam} \\
  --primers {primers} \\
  --output {output_bam} \\
  {hard_clip} \\
  {other_options}""".format(
            tmp_dir=config.param('fgbio_trim_primers', 'tmp_dir'),
            java_other_options=config.param('fgbio_trim_primers', 'java_other_options'),
            ram=config.param('fgbio_trim_primers', 'ram'),
            input_bam=input_bam,
            primers=config.param('fgbio_trim_primers', 'primers'),
            output_bam=output_bam,
            hard_clip="-H" if hard_clip else "",
            other_options=config.param('fgbio_trim_primers', 'other_options')
            )
        )

def groupreadsbyumi(
    input_bam,
    output=None,
    metrics_dir=None
    ):

    inputs = [input_bam]
    outputs = [output]
    return Job(
        inputs,
        outputs,
        [
            ['fgbio_groupreadsbyumi', 'module_java'],
            ['fgbio_groupreadsbyumi', 'module_fgbio']
        ],

        command="""\
java -Djava.io.tmpdir={tmp_dir} {java_other_options} -Xmx{ram} -jar $FGBIO_JAR {fgbio_options}GroupReadsByUmi \\
  -i {input_bam}{output_file}{metrics_dir} \\
  {other_options}""".format(
        tmp_dir=config.param('fgbio_groupreadsbyumi', 'tmp_dir'),
        java_other_options=config.param('fgbio_groupreadsbyumi', 'java_other_options'),
        ram=config.param('fgbio_groupreadsbyumi', 'ram'),
        input_bam=input_bam,
        output_file=" -o " + output if output else "",
        metrics_dir=" -f " + metrics_dir if metrics_dir else "",
        other_options=config.param('fgbio_groupreadsbyumi', 'other_options') if config.param('fgbio_groupreadsbyumi', 'other_options',required=False) else "",
        fgbio_options=config.param('fgbio_groupreadsbyumi', 'fgbio_options') + " " if config.param('fgbio_groupreadsbyumi', 'fgbio_options',required=False) else " "
        )
    )

def callduplexconsensusreads(
    input,
    output
    ):

    inputs = [input]
    outputs = [output]
    return Job(
        inputs,
        outputs,
        [
            ['fgbio_callduplexconsensusreads', 'module_java'],
            ['fgbio_callduplexconsensusreads', 'module_fgbio']
        ],

        command="""\
java -Djava.io.tmpdir={tmp_dir} {java_other_options} -Xmx{ram} -jar $FGBIO_JAR GroupReadsByUmi \\
  -i {input}{output_file} \\
  {other_options}""".format(
        tmp_dir=config.param('fgbio_callduplexconsensusreads', 'tmp_dir'),
        java_other_options=config.param('fgbio_callduplexconsensusreads', 'java_other_options'),
        ram=config.param('fgbio_callduplexconsensusreads', 'ram'),
        input=input,
        output_file=" -o " + output if output else "",
        other_options=config.param('fgbio_callduplexconsensusreads', 'other_options') if config.param('fgbio_callduplexconsensusreads', 'other_options',required=False) else ""
        )
    )

def callmolecularconsensusreads(
    input,
    output,
    group_name=None
    ):

    inputs = [input]
    outputs = [output]
    return Job(
        inputs,
        outputs,
        [
            ['fgbio_callmolecularconsensusreads', 'module_java'],
            ['fgbio_callmolecularconsensusreads', 'module_fgbio']
        ],

        command="""\
java -Djava.io.tmpdir={tmp_dir} {java_other_options} -Xmx{ram} -jar $FGBIO_JAR CallMolecularConsensusReads \\
-i {input}{output_file}{group_name} \\
  {other_options}{threads}""".format(
        tmp_dir=config.param('fgbio_callmolecularconsensusreads', 'tmp_dir'),
        java_other_options=config.param('fgbio_callmolecularconsensusreads', 'java_other_options'),
        ram=config.param('fgbio_callmolecularconsensusreads', 'ram'),
        input=input,
        output_file=" -o " + output if output else "",
        group_name=" -R " + group_name if group_name else "",
	threads =" --threads " +  config.param('fgbio_callmolecularconsensusreads', 'threads') if config.param('fgbio_callmolecularconsensusreads', 'threads') else "",
        other_options=config.param('fgbio_callmolecularconsensusreads', 'other_options') if config.param('fgbio_callmolecularconsensusreads', 'other_options',required=False) else ""
        )
    )

def filterconsensusreads(
    input,
    output
    ):

    inputs = [input]
    outputs = [output]
    return Job(
        inputs,
        outputs,
        [
            ['fgbio_filterconsensusreads', 'module_java'],
            ['fgbio_filterconsensusreads', 'module_fgbio']
        ],
        command="""\
java -Djava.io.tmpdir={tmp_dir} {java_other_options} -Xmx{ram} -jar $FGBIO_JAR FilterConsensusReads \\
  -i {input}{output_file}{reference} \\
  {other_options}""".format(
        tmp_dir=config.param('fgbio_filterconsensusreads', 'tmp_dir'),
        java_other_options=config.param('fgbio_filterconsensusreads', 'java_other_options'),
        ram=config.param('fgbio_filterconsensusreads', 'ram'),
        input=input,
        output_file=" -o " + output if output else "",
        reference=" -r " + config.param('fgbio_filterconsensusreads', 'genome_fasta') if config.param('fgbio_filterconsensusreads', 'genome_fasta',required=False) else "",
        other_options=config.param('fgbio_filterconsensusreads', 'other_options') if config.param('fgbio_filterconsensusreads', 'other_options',required=False) else ""
        )
    )

def extractumisfrombam(
    input_bam,
    output_bam
    ):

    inputs = [input_bam]
    outputs = [output_bam]
    return Job(
        inputs,
        outputs,
        [
            ['fgbio_extractumisfrombam', 'module_java'],
            ['fgbio_extractumisfrombam', 'module_fgbio']
        ],

        command="""\
java -Djava.io.tmpdir={tmp_dir} {java_other_options} -Xmx{ram} -jar $FGBIO_JAR ExtractUmisFromBam \\
  --input {input_bam} \\
  --output {output_bam} \\
  {other_options}""".format(
        tmp_dir=config.param('fgbio_extractumisfrombam', 'tmp_dir'),
        java_other_options=config.param('fgbio_extractumisfrombam', 'java_other_options'),
        ram=config.param('fgbio_extractumisfrombam', 'ram'),
        input_bam=input_bam,
        output_bam=output_bam,
	other_options=config.param('fgbio_extractumisfrombam', 'other_options') if config.param('fgbio_extractumisfrombam', 'other_options',required=False) else ""
        ),
        removable_files=[output_bam]
    )

def sortbam(
    input,
    output
    ):

    inputs = [input]
    outputs = [output]
    return Job(
        inputs,
        outputs,
        [
            ['fgbio_callmolecularconsensusreads', 'module_java'],
            ['fgbio_callmolecularconsensusreads', 'module_fgbio']
        ],

        command="""\
java -Djava.io.tmpdir={tmp_dir} {java_other_options} -Xmx{ram} -jar $FGBIO_JAR SortBam \\
-i {input}{output_file} \\
  {other_options}""".format(
        tmp_dir=config.param('fgbio_sort', 'tmp_dir'),
        java_other_options=config.param('fgbio_sort', 'java_other_options'),
        ram=config.param('fgbio_sort', 'ram'),
        input=input if input else "/dev/stdin",
        output_file=" -o " + output if output else "/dev/stdout",
        other_options=config.param('fgbio_sort', 'other_options') if config.param('fgbio_sort', 'other_options',required=False) else ""
        )
    )

