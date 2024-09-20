# CHIP COLCOT scripts

The present repository contains source code and guidelines for execution of implemented scripts that were performed in Digital Research Alliance of Canada's Narval cluster. Sequencing data analysis was performed with a workflow implemented in GenPipes as per Paragon’s instructions. 

# Table of Contents

-   [CHIP COLCOT
    scripts](#chip-colcot-scripts)
    -   [Clonal Hematopoiesis Sequencing Data
        Processing](#clonal-hematopoiesis-sequencing-data-processing)
        -   [Clonal Hematopoiesis Genpipes
            Pipeline](#clonal-hematopoiesis-genpipes-pipeline)
    -   [Variant calls and filtering
        workflow](#variant-calls-and-filtering-workflow)
        -   [1 - Vardict](#vardict)
        -   [2 - Variant decompose and normalize using
            vt](#variant-decompose-and-normalize-using-vt)
        -   [3 - Variants annotations VEP and
            Annovar](#variants-annotations-vep-and-annovar)
        -   [4 - Chip filtering](#chip-filtering)
        -   [5 - Process
            variants](#process-variants)
    -   [Report and output
        generation](#report-and-output-generation)
        -   [1 - Binomial tests
            annotations](#binomial-tests-annotations)
        -   [2 - Variants
            filters](#variants-filters)
        -   [3 - CHIP hotspots
            filtering](#chip-hotspots-filtering)
    -   [Depth of coverage
        extraction](#depth-of-coverage-extraction)
        -   [1 - Vardict in pileup
            mode](#vardict-in-pileup-mode)
        -   [2 - Variant decompose and normalize using
            vt](#variant-decompose-and-normalize-using-vt-1)
        -   [3 - Index variants using
            tabix](#index-variants-using-tabix)
        -   [4 - Compile
            results](#compile-results)
        -   [5 - Generate output in wide
            format.](#generate-output-in-wide-format.)


## Clonal Hematopoiesis Sequencing Data Processing

### Clonal Hematopoiesis Genpipes[1] Pipeline

We implemented the recommended Paragon Unique Molecular Identifiers (UMI) 12 step analysis. After demultiplexing sequencing data, trailing sequencing adapters in the R1 and R2 reads are trimmed using cutadapt v.4.6. Trimmed fastq reads are then converted to unmapped alignment (BAM) files. The 16bp UMI is removed from BAM files to be stored as an alignment tag. Template reads are reconverted to fastq format and mapped against the human genome reference GRCh37 using bwa-mem v.0.7.17. The mapped and unmapped tagged reads are merged and the alignment file (with suffix mapped.umi.bam) is created. The resulting alignments are processed using Fulcrum Genomics’ consensus calling workflow provided by the toolset fgbio v2.2.1. Aligned reads are grouped as read families that appeared to have come from the same original molecule by adjacency strategy and allowing one edit between UMIs (only bases with mapping quality >= 20 are kept). Fgbio CallMolecularConsensusReads is then used to create consensus reads by combining the evidence from reads with the same unique molecular tag: at least 3 reads are required to produce a consensus base. The consensus reads are mapped against the human genome reference using bwa-mem v.0.7.17. Alignment files are filtered to exclude alignments with insert sizes <75 and > 305.


#### Example run

``` bash
python GenPipes/pipelines/clonal_hematopoiesis/clonal_hematopoiesis.py \
-r readset.tsv \
-c GenPipes/pipelines/clonal_hematopoiesis/clonal_hematopoiesis.base.ini \
GenPipes/pipelines/clonal_hematopoiesis/clonal_hematopoiesis.beluga.ini \
GenPipes/pipelines/clonal_hematopoiesis/clonal_hematopoiesis.novaseq.ini \
-s 2-9,12,15 \
--no-json > toRun.sh

bash toRun.sh
```

## Variant calls and filtering workflow


### 1 - Vardict

The vardict v1.8.2 variant caller is used to obtain variant calls.

``` bash
cd ./vardict/
bash execute_vardict.sh
```

### 2 - Variant decompose and normalize using vt

The vt tool v.0.57 is used to decompose and normalize complex variants in VCF files.

``` bash
cd ./vt_decompose_normalize
./execute_vt_decompose_normalize.sh
```


### 3 - Variants annotations VEP and Annovar

Variant consequence prediction was performed using Annovar v.2020-06-08.  Variants observed in Topmed and known sequencing artifacts were flagged using vep v.110.1

``` bash
cd ./vep_annovar/
./execute_vep_annovar.sh
```


### 4 - Chip filtering

Variants were flagged as CHIP (or whitelisted) if they match a pre-specified list of putative CHIP variants as in Vlasschaert et al.
R version 4.1.2 was used and required packages are listed in [R_session_info.txt](https://github.com/pgxcentre/COLCOT-CH-scripts/blob/main/R_session_info.txt)

``` bash
cd ./chip_filtering
./execute_chip_filtering.sh

```

### 5 - Process variants

Execute a python script was used to perform modifications including:

- Compute Variant allele fraction (VAF) using alternate allele (AD) divided by total read depth (DP).

- Remove variants with a total read depth (DP) of < 20.

- Remove known sequencing artifacts reported by Busque et al. and flagged by vep. (TP53:NM_000546.6:c.1129A>C, TP53:NM_000546.6:c.215C>G, ASXL1:NM_015338.6:c.708G>C , ASXL1:NM_015338.6:c.2444T>C). For NM_015338.6:c.1934dupG variants are removed if detected allele fraction is less than 6%.
- Remove frameshift mutations if they occurred in homo-polymer repeats (5 consecutive reads of the same nucleotide) unless there was a total 10 or more supporting reads and a VAF>8% for these variants.
- Perform post-modifications to the whitelist filters as follows:

  –	Exclude variants from ASXL1 which are not located in exon 11 and 12. These should have wl.exception set to FALSE, but whitelist set to TRUE.

  –	Exclude synonymous variants that were marked as whitelist, because of the way the scripts looks for stop gained / loss (loss of function). It looks with a regular expression for a X in the NonsynOI column. But if there is a synonymous variant which changes a STOP for a STOP (for example, X123X means stop à 123 for a stop).

  – Variants originally marked for manual verification are marked as whitelist.

``` bash
cd ./process_variants
./process_vlasschaert_variants.sh
```

## Report and output generation

### 1 - Binomial tests annotations

Binomial tests annotations added to show if a variant significantly deviates from the expected distribution for a germline allele (defined as a p-value from a binomial test of less than 0.01 assuming a probability of success in a single Bernoulli experiment of 0.5 and using the alternate allele read count (VD) as the number of successes and the total depth (DP) as the number of trials).

``` bash
cd report

Rscript.sh binomial_tests.Rscript \
--annovar_path ../process_variants/all_samples.varsOI.allvariants.processed.tsv \
--output_prefix all_samples.varsOI.allvariants.processed

```

### 2 - Variants filters

Therefore, an R script available at Gitlab was used to filter putative CHIP variants marked as whitelist with a total read depth (DP) >= 20.

``` bash
./Rscript.sh filter_variants_DP.Rscript \
--annovar_path all_samples.varsOI.allvariants.processed_annotated_binomial.tsv \
--data_keys NovaSeq_Samples_2_filtered.csv \
--output_prefix all_samples.varsOI.allvariants.processed.filtered

```


### 3 - CHIP hotspots filtering

Variants observed in more than 10 participants were assessed for associations with age and a common genetic variant in the TERT promoter (rs7705526) using Firth’s bias reduction logistic regression provided by the R logistf package. Missense CHIP putative variants listed in Jaiswal et al.2 or chip hotspots associated to age or TERT in the UK Biobank2 were excluded from association tests.

``` bash
cd associations_age_tert

./Rscript.sh associations_all_variants_processed.Rscript \
-a ../all_samples.varsOI.allvariants.processed.filtered_colcot_filtered_DP_20.tsv \
-g ./rs7705526_imputed.raw -p ./COLCOT_age.csv \
-s ./chip_samples_baseline.txt \
-o ./all_samples.varsOI.allvariants.processed.filtered_DP_20_firth_w_exclusions \
-f -m ./ressources/whitelist_filter_files/CHIP_ssociation_exceptions_2024-07-05.txt

```

Prepare file with minimal information for descriptive statistics and tests.

``` bash

./Rscript.sh prepare_final_chip_file.Rscript \
--annovar_path ./all_samples.varsOI.allvariants.processed.filtered_DP_20_firth_w_exclusions_w_association_exclusions.tsv \
--old_names 'PGx.Sample.ID,VAF,variant_Id' \
--new_names 'PGx_Sample_ID,vaf,variant_ID' \
--output_prefix colcot_CHIP_filtered_DP_20_whitelist \
--select_variants 'PGx_Sample_ID,SAMPLE,SYMBOL,REF,ALT,DP,VD,vaf,AF,variant_ID'

```

## Depth of coverage extraction

Extract coverage and variant depth for CHIP variants identified within 848 COLCOT patients with DNA available and sequenced at baseline and at the end of the study, including variants not detected or filtered due to allelic fraction threshold in vardict (AF<0.001). 

### 1 - Vardict in pileup mode 

The vardict v1.8.2 variant caller is used to obtain records in vcf format of all possible variants in mpileup mode.

``` bash
cd ./vardict_all/
bash execute_vardict.sh
```

### 2 - Variant decompose and normalize using vt

The vt tool v.0.57 is used to decompose and normalize complex variants in VCF files.

``` bash
cd extract_dp_all_variants

find ../vardict_all/output/ \
-name "*AF_consensus.vardict.vcf.gz" > all_vardict_to_process

parallel --results vt_normalize \
--joblog vt_normalize.log -n 1 \
-j 20 ./vt_decompose_normalize.sh \
${source_path}/{} ${ref_fasta} output/{} :::: all_vardict_to_process
```

### 3 - Index variants using tabix

Extract variants detected at every known CHIP position.

``` bash
parallel --rpl '{%(.+?)} s/$$1$//;' \
--results tabix_counts \
--joblog tabix.log -n 1 -j 12 ./tabix.sh \
  {} ./unique_variant_positions.tsv {%vcf.gz.vcf.gz}tsv ::: output/*AF_consensus.vardict.vcf.gz.vcf.gz
```

### 4 - Compile results 

Add pileup counts to original variants file, generate a file with minimum information for variants called and not called in any visit.

``` bash
./Rscript.sh add_counts_variants_per_visit.Rscript \
--annovar_path ../associations_age_tert/colcot_CHIP_filtered_DP_20_whitelist.tsv \
--data_keys ./NovaSeq_Samples_2_filtered.csv \
--counts_path output --pattern ".AF_consensus.vardict.tsv$" \
--output_prefix colcot_CHIP_filtered_DP_20_whitelist
```

### 5 - Generate output in wide format.

``` bash
./Rscript.sh convert_to_wide.Rscript \
--annovar_path colcot_CHIP_filtered_DP_20_whitelist_both_visits_w_counts.tsv \
--output_prefix SGR-2910 > _convert_to_wide.log 2>&1
```


## References

[1] Bourgey M, Dali R, Eveleigh R, Chen KC, Letourneau L, Fillon J, Michaud M, Caron M, Sandoval J, Lefebvre F, et al. GenPipes: an open-source framework for distributed and scalable genomic analyses. GigaScience. 2019;8. doi: 10.1093/gigascience/giz037