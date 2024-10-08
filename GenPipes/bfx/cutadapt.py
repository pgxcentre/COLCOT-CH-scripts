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

# MUGQIC Modules
from core.config import *
from core.job import *

log = logging.getLogger(__name__)

def trim(input1, input2, prefix, adapter_3p_fwd, adapter_3p_rev, metrics_dir=""):
    # output_pair1 = prefix + ".trim.pair1.fastq.gz"
    # output_pair2 = prefix + ".trim.pair2.fastq.gz"

    if input2:  # Paired end reads
        output_pair1 = prefix + ".trim.pair1.fastq.gz"
        output_pair2 = prefix + ".trim.pair2.fastq.gz"
        inputs = [input1, input2]
        output = [output_pair1, output_pair2]
    else:   # Single end reads
        output_pair1 = prefix + ".trim.single.fastq.gz"
        output_pair2 = None
        inputs = [input1]
        output = [output_pair1]

    return Job(
        inputs,
        output,
        [
            ['cutadapt', 'module_cutadapt']
        ],

        command="""\
cutadapt {adapter_5p_fwd} \\
  {adapter_5p_rev} \\
  {adapter_3p_fwd} \\
  {adapter_3p_rev} \\
  {nthread} \\
  {options} \\
  {output_fwd} \\
  {output_rev} \\
  {inputs}{log}""".format(
      adapter_5p_fwd="-g " + config.param('cutadapt', 'adapter_5p_fwd') if config.param('cutadapt', 'adapter_5p_fwd') else "",
      adapter_5p_rev="-G " + config.param('cutadapt', 'adapter_5p_rev') if input2 and config.param('cutadapt', 'adapter_5p_fwd') else "",
      adapter_3p_fwd="-a " + adapter_3p_fwd,
      adapter_3p_rev="-A " + adapter_3p_rev if input2 else "",
      options=config.param('cutadapt', 'options'),
      nthread="-j " + config.param('cutadapt', 'threads'),
      inputs=" \\\n  ".join(inputs),
      output_fwd="-o " + output_pair1,
      output_rev="-p " + output_pair2 if input2 else "",
      log=" > " + os.path.join(metrics_dir, re.sub(r"(\w.*?)\.(\w.*?)$",r"\1.log", os.path.basename(output_pair1))) if metrics_dir !="" else "" 
      ),
    )
