# Variant calls filtering

Post-variant call filtering was implemented in the script whitelist_filter_rscript.R described by Vlasschaert et al[1]. (https://github.com/briansha/Annovar_Whitelist_Filter_WDL) with required resources extracted from Terra Workspaces (https://app.terra.bio/#workspaces/terra-outreach/CHIP-Detection-Mutect2), including : NEJM_2017_genes_01262020.txt, CHIP_splice_vars_agb_01262020.txt, CHIP_nonsense_FS_vars_agb_01262020.txt, and CHIP_missense_vars_agb_01262020.txt.

The script was cloned and modified to account for minor script errors, including the detection of nucleotide duplication introducing an immediate stop codon without a frameshift (denoted by “*”). 


## Referencez


[1] Vlasschaert C, Mack T, Heimlich JB, Niroula A, Uddin MM, Weinstock J, Sharber B, Silver AJ, Xu Y, Savona M, et al. A practical approach to curate clonal hematopoiesis of indeterminate potential in human genetic data sets. Blood. 2023;141:2214-2223. doi: 10.1182/blood.2022018825