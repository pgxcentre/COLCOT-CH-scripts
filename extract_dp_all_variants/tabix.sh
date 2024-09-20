#!/usr/bin/env bash

set -eo pipefail
module purge
module load mhipgx/htslib/1.19

main() {
    local input=$1
    local target=$2
    local output=$3

if [[ -z "$output" ]]; then
  output=${input/.vcf/.target.vcf}
fi;

tabix ${input} -T ${target} > ${output}

}

main "$@"
