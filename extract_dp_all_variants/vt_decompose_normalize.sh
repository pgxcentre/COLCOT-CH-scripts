#!/usr/bin/env bash

set -eo pipefail

module purge
module load mhipgx/htslib/1.19 mugqic/vt/0.57

main() {
local input_vcf=$1
local assembly_fasta=$2
local output=$3

cat ${input_vcf} | vt decompose -s - | vt normalize -r ${assembly_fasta} - | \
bgzip -cf \
 > \
${output}.vcf.gz && tabix -pvcf ${output}.vcf.gz

}

main "$@"
