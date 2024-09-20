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

# MUGQIC Modules
from posixpath import join
from core.config import *
from core.job import *

#vep  --cache --offline --force_overwrite -i $fi --dir_cache $VEP_CACHE --dir_plugins $VEP_PLUGINS --fasta $fasta --assembly GRCh37 -o vep/$sname.vep_out --stats_file vep/$sname.vep_out.stats.html --force_overwrite --fork 8 --format vcf 
#--everything --check_existing --clin_sig_allele 0 -distance 5000 --hgvs --mane --regulatory --symbol --transcript_version --tsl --var_synonyms --plugin FATHMM_MKL,$VEP_CACHE/fathmm-MKL/fathmm-MKL_Current.tab.gz --custom $VEP_CACHE/clinvar/clinvar.vcf.gz,ClinVar,vcf,exact,0,CLNSIG,CLNREVSTAT,CLNDN
def annotate(input, output=None):
    plugins=""
    if config.param('ensembl_vep', 'plugins', required=False):
        plugins_list=config.param('ensembl_vep', 'plugins').split(":")
        plugins='{1}{0}'.format(" --plugin ".join(plugins_list), " --plugin ")
    custom=""
    if config.param('ensembl_vep', 'custom', required=False):
        custom_list=config.param('ensembl_vep', 'custom').split(":")
        custom='{1}{0}'.format(" --custom ".join(custom_list), " --custom ")
    if config.param('ensembl_vep', 'output_format', required=False):
	output_format =  config.param('ensembl_vep', 'output_format', required=False)
    else:
	output_format = 'vcf'
    # Job command
    return Job(
        [input],
        [output],
        [
        ['ensembl_vep', 'module_vep'],
        ],
        command="""\
vep \\
--cache --offline --force_overwrite \\
-i {input} \\
--dir_cache $VEP_CACHE \\
--dir_plugins $VEP_PLUGINS \\
--fasta {reference_fasta} \\
--assembly {assembly} \\
-o {output} \\
--{output_format} \\
--stats_file {output}_summary.html \\
--fork {threads} \\
--format vcf {other_options}{plugins}{custom}""".format(
        input=input,
        assembly=config.param('ensembl_vep', 'assembly'),
        reference_fasta=config.param('ensembl_vep', 'genome_fasta', type='filepath'),
        other_options=config.param('ensembl_vep', 'other_options') if config.param('ensembl_vep', 'other_options', required=False) else "--everything",
        threads=config.param('ensembl_vep', 'threads', required=False) if config.param('ensembl_vep', 'threads', required=False) else 1,
        plugins=plugins,
        custom=custom,
        output=output,
	output_format=output_format
        )
    )

def filter(input, output=None):

    filter_list = config.param('ensembl_vep', 'filters').split(":")
    filters = '{1}{0}'.format(" -filter ".join(filter_list), " -filter ")
    # Job command
    return Job(
        [input],
        [output],
        [
        ['ensembl_vep', 'module_vep'],
        ],
        command="""\
filter_vep --force_overwrite -i {input} \\
-o {output}{filters}""".format(
        input=input,
        filters=filters,
        output=output
        )
    )
