#!/usr/bin/env bash

# STEP 3


set -eo pipefail

module load apptainer/1.1.8 \
            mhipgx/htslib/1.9 \
            mhipgx/annovar/2021_12_09

execute_vep_annovar() {
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

    (
        apptainer run -C -B /lustre07/scratch/"$USER" -B /lustre07/scratch/"$USER":/scratch/"$USER" -W "$PWD" ~/apptainer_sif/ensembl-vep-110.1.sif vep \
            --cache \
            --offline \
            --force_overwrite \
            -i "$PWD"/"$vcf_file" \
            --dir_cache "$PWD"/vep_cache \
            --dir_plugins "$PWD"/vep_plugins \
            --fasta "$PWD"/ressources/Homo_sapiens.GRCh37.fa \
            --assembly GRCh37 \
            --refseq \
            -o STDOUT \
            --vcf \
            --fork 1 \
            --format vcf \
            --CACHE_VERSION 110 \
            --custom "$PWD"/ressources/InheritedCausesofClonalHematopoiesisin97691TOPMedWholeGenomes.GRCh37.vcf.gz,TopmedClonal,vcf,exact,0 \
            --custom "$PWD"/ressources/Busque_artifacts_GRCh37.vcf.gz,Busque_possible_artifact,vcf,exact \
        | bgzip -cf \
        > "$out_dir"/"$out_prefix".vep.vcf.gz
    ) 2> "$out_dir"/"$out_prefix".vep.log

    # Indexing
    tabix -pvcf "$out_dir"/"$out_prefix".vep.vcf.gz

    # Annotating with Annovar
    (
        table_annovar.pl "$out_dir"/"$out_prefix".vep.vcf.gz "$HUMANDB" \
            -buildver hg19 \
            -out "$out_dir"/"$out_prefix" \
            -remove \
            -protocol refGene,cosmic70 \
            -operation g,f \
            -nastring . \
            --convertarg '-include' \
            -vcfinput
    ) > "$out_dir"/"$out_prefix".annovar.log 2>&1
}
export -f execute_vep_annovar

main() {
    if [[ $# -eq 0 ]]; then
        # Finding all the VCF files
        find ../vt_decompose_normalize -type f -name "*.vcf.gz" -print0 \
        | shuf -z \
        | split -t '\0' -l 278 --additional-suffix .txt - ./vep_annovar_input_

        # Executing by sumbitting the script
        for fn in ./vep_annovar_input_*.txt; do
            ./execute_vep_annovar.sh "$fn"
        done
    else
        for fn in "$@"; do
            local suffix
            suffix=$(echo "${fn%\.txt}" | grep -o -E "_[a-z]+$")

            parallel --null --joblog ./parallel"$suffix".log --resume-failed -n1 -j30 \
                execute_vep_annovar \
                :::: "$fn"
        done
    fi
}

main "$@"
