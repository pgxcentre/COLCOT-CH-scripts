#!/usr/bin/env bash

# STEP 1
# ./execute_vardict.sh

set -eo pipefail

export MUGQIC_INSTALL_HOME=/cvmfs/soft.mugqic/CentOS6
module use "$MUGQIC_INSTALL_HOME"/modulefiles
module load mugqic/java/openjdk-jdk1.8.0_72 \
            mhipgx/VarDictJava/1.8.3 \
            samtools/1.16.1 \
            mugqic/perl/5.22.1 \
            mhipgx/R_Bioconductor/4.1.2_3.14 \
            mhipgx/htslib/1.9

execute_vardict() {
    local bam_file=$1

    # Finding the output dir
    local out_dir=./output

    # Finding the sample ID
    local sample_id
    sample_id=$(basename "${bam_file%\.sorted\.filtered\.bam}")

    (
        VarDict \
            -G /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh37/genome/Homo_sapiens.GRCh37.fa \
            -N "$sample_id" \
            -b "$bam_file" \
            -a 10:0.95 -f 0.001 -q 20 -r 1 -k 1 -P 2 -I 50 -c 1 -S 2 -E 3 -g 4 -U -L 51 -p \
            PGD628_v3.ampInfo.bed \
        | "$VARDICT_BIN"/teststrandbias.R \
        | "$VARDICT_BIN"/var2vcf_valid.pl \
            -A \
            -N "$sample_id" \
            -E \
            -f 0.001 \
        | bgzip -c \
        > "$out_dir"/"$sample_id".AF_consensus.vardict.vcf.gz
    ) 2> "$out_dir"/"$sample_id".AF_consensus.vardict.log

    # Indexing
    tabix -pvcf "$out_dir"/"$sample_id".AF_consensus.vardict.vcf.gz
}
export -f execute_vardict

main() {
    if [[ $# -eq 0 ]]; then
        # Finding all the BAM files
        find -L ../bam -type f -name "*.bam" -print0 \
        | shuf -z \
        | split -t '\0' -l 278 --additional-suffix .txt - ./vardict_input_

        # Executing VarDict by submitting the script
        for fn in ./vardict_input_*.txt; do
            ./execute_vardict.sh "$fn"
        done
    else
        for fn in "$@"; do
            local suffix
            suffix=$(echo "${fn%\.txt}" | grep -o -E "_[a-z]+$")

            parallel --null --joblog ./parallel"$suffix".log --resume-failed -n1 -j30 \
                execute_vardict \
                :::: "$fn"
        done
    fi
}

main "$@"
