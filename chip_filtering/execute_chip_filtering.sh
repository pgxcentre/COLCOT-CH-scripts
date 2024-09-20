#!/usr/bin/env bash

# STEP 4

set -eo pipefail

module load mhipgx/R_Bioconductor/4.1.2_3.14

execute_chip_filtering() {
    local annovar_file=$1

    # Finding the output dir
    local out_dir=./output

    # The prefix
    local sample_id
    sample_id=$(basename "${annovar_file%\.txt}")

    Rscript --vanilla ./whitelist_filter_vlasschaert_rscript.R \
        ./ressources/whitelist_filter_files \
        "$annovar_file" \
        "$out_dir" \
        "$sample_id"
}
export -f execute_chip_filtering

main() {
    # Finding all the annovar files
    find ../vep_annovar -type f -name "*.hg19_multianno.txt" -print0 | sort -z > ./chip_filtering_input.txt

    # Executing annotation
    parallel --null --joblog ./parallel.log --resume-failed -n1 -j30 \
        execute_chip_filtering \
        :::: ./chip_filtering_input.txt
}

main "$@"
