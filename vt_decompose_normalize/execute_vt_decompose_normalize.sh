#!/usr/bin/env bash

# STEP 2

set -eo pipefail

export MUGQIC_INSTALL_HOME=/cvmfs/soft.mugqic/CentOS6
module use "$MUGQIC_INSTALL_HOME"/modulefiles

module load mhipgx/htslib/1.9 \
            mugqic/vt/0.57

execute_vt_decompose_normalize() {
    local vcf_file=$1

    # First, we check if we have variants
    local nb_variants
    nb_variants=$(zcat "$vcf_file" | grep -c -v '^#')
    if [[ $nb_variants -eq 0 ]]; then
        echo "${vcf_file}: no variants" >&2
        return 1
    fi

    # Finding the output dir
    local out_dir=./output

    # The prefix
    local out_prefix
    out_prefix=$(basename "${vcf_file%\.vcf\.gz}")

    # Executing
    (
        zcat "$vcf_file" \
        | vt decompose -s - \
        | vt normalize -r /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh37/genome/Homo_sapiens.GRCh37.fa - \
        | bgzip -cf \
        > "$out_dir"/"$out_prefix".vt.vcf.gz
    ) > "$out_dir"/"$out_prefix".vt.log 2>&1

    # Indexing
    tabix -pvcf "$out_dir"/"$out_prefix".vt.vcf.gz
}
export -f execute_vt_decompose_normalize

main() {
    # Finding all the VCF files
    find ../vardict -type f -name "*.vcf.gz" -print0 | sort -z > ./vt_decompose_normalize_input.txt

    # Executing annotation
    parallel --null --joblog ./parallel.log --resume-failed -n1 -j30 \
        execute_vt_decompose_normalize \
        :::: ./vt_decompose_normalize_input.txt
}

main "$@"
