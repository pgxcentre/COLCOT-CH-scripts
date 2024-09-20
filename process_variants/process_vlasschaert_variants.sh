#!/usr/bin/env bash

# STEP 5

set -eo pipefail

ml python/3.11.2

compare_header() {
    local fn=$1

    diff -q ./header.csv <(head -n1 "$fn")
}
export -f compare_header

main() {
    # Finding all the files
    find ../chip_filtering -type f -name "*.allvariants.csv" -print0 | sort -z > ./process_variants_input.txt

    # The first file of the list
    local first_file
    first_file=$(head -z -n1 ./process_variants_input.txt | sed -e 's/\x0//')

    # The header of the first file
    head -n1 "$first_file" > ./header.csv

    # Checking the header
    parallel --null --joblog ./parallel.log --resume-failed -n1 -j30 \
        compare_header \
        :::: ./process_variants_input.txt

    # Concatenating the files (adding the input filename)
    (
        paste -d, <(head -n1 "$first_file") <(echo '"Input.File"')

        while IFS= read -r -d $'\0' fn; do
            sed -e 1d -e 's+$+,"'"$fn"'"+' "$fn"
        done < ./process_variants_input.txt
    ) > ./all_samples.varsOI.allvariants.csv

    # Processing the variants
    ./process_vlasschaert_variants.py \
        --vaf 0.0 \
        --vd 1 \
        --jak2-ok-manual F537_K539delinsL \
        --vep-header "Allele|Consequence|IMPACT|SYMBOL|Gene|Feature_type|Feature|BIOTYPE|EXON|INTRON|HGVSc|HGVSp|cDNA_position|CDS_position|Protein_position|Amino_acids|Codons|Existing_variation|DISTANCE|STRAND|FLAGS|SYMBOL_SOURCE|HGNC_ID|REFSEQ_MATCH|REFSEQ_OFFSET|GIVEN_REF|USED_REF|BAM_EDIT|SOURCE|TopmedClonal|Busque_possible_artifact" \
        --reference /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh37/genome/Homo_sapiens.GRCh37.fa \
        ./all_samples.varsOI.allvariants.csv \
        ./all_samples.varsOI.allvariants.processed.tsv \
    2> ./process_vlasschaert_variants.log

    rm -f ./header.csv
}

main "$@"
