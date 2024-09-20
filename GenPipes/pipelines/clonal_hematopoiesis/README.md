[TOC]


Clonal Hematopoiesis Amplicon Pipeline
=================
The Clonal Hematopoiesis Pipeline is based on the DNA-Seq pipeline and identifies somatic variants
within target sequencing data. Somatic variants are identified by allowing for varying ploidy in
the form of allele fractions for each variant. If Unique Molecular Identifiers are detected the
fgbio CallMolecularConsensusReads tool is used. Recalibration is not run either because typically,
these datasets are targeted with amplicons or custom capture which render recalibration useless.

Data Analysis Recommendations for Illumina Paired-end Sequencing using CleanPlex UMI Panels follows:

1- Trim trailing adapters in the R1 and R2 reads using cutadapt
    cutadapt -g CCTACACGACGCTCTTCCGATCT \
    -a AGATCGGAAGAGCACACGTCTGAA \
    -A AGATCGGAAGAGCGTCGTGTAGG \
    -G TTCAGACGTGTGCTCTTCCGATCT -o SampleID.1.fastq -p SampleID.2.fastq SampleID.rawreads.1.fastq SampleID.rawreads.2.fastq -e 0.1 -O 9 -m 20 -n 2
    Only when the panel is GCP53, cutadapt needs to be run for a second time:
    cutadapt -g GTTATTGGACGTTTA \
    -a TTACCTTTTATGAAG \
    -A TAAACGTCCAATAAC \
    -G CTTCATAAAAGGTAA \
    -e 0.1 -O 9 -m 20 -n 2 \
    -o SampleID.trimmed.1.fastq -p SampleID.trimmed.1.fastq SampleID.1.fastq  SampleID.2.fastq

2- Step 1: Convert the trimmed fastq reads to unmapped BAM
    java -jar picard.jar FastqToSam F1=SampleID.trimmed.1.fastq F2=SampleID.trimmed.2.fastq         SM=SampleID O=SampleID.step1-unmapped.bam

3- Step 2: Extract the UMI information to store as a tag
    java -jar fgbio.jar ExtractUmisFromBam --input=SampleID.step1-unmapped.bam --output=SampleID.step2-unmapped-umi.bam --read-structure=16M+T 16M+T --single-tag=RX --molecular-index-tags=ZA ZB

4- Step 3: Convert the template read to fastq
    java -jar picard.jar SamToFastq I=SampleID.step2-unmapped-umi.bam F=SampleID.step3-R1Template.fastq F2=SampleID.step3-R2Template.fastq

5- Step 4: Map the template reads against the human genome
    bwa mem -t No_threads_computation  -M hg19.fa SampleID.step3-R1Template.fastq SampleID.step3-R2Template.fastq > SampleID.step4-MappedTemplate.bam

6- Step 5: Merge the mapped reads with the UMI info
    java -jar picard.jar MergeBamAlignment UNMAPPED=SampleID.step2-unmapped-umi.bam ALIGNED=SampleID.step4-MappedTemplate.bam R=hg19.fa O=SampleID.step5-Mapped-UMI.bam

7- Step 6: Group reads by UMI
    java -jar fgbio.jar GroupReadsByUmi --input=SampleID.step5-Mapped-UMI.bam --output=SampleID.step6-Mapped-GroupedByUMI.bam --strategy=adjacency --edits=1 --min-map=10

8- Step 7: Call consensus on grouped reads
    java -jar fgbio.jar CallMolecularConsensusReads --input=SampleID.step6-Mapped-GroupedByUMI.bam --output=SampleID.step7-ConsensusReads-Unmapped.bam --min-reads=3

9- Step 8: Convert the consensus reads to fastq
    java -jar picard.jar SamToFastq I=SampleID.step7-ConsensusReads-Unmapped.bam F=SampleID.step8-R1Consensus.fastq F2=SampleID.step8-R2Consensus.fastq

10- Step 9: Map the consensus reads against the human genome
    bwa mem -t No_threads_computation -M hg19.fa SampleID.step8-R1Consensus.fastq SampleID.step8-R2Consensus.fastq > SampleID.step9-Mapped-Consensus.bam

11- Step 10: Merge the unmapped and the mapped consensus files
    java -jar picard.jar MergeBamAlignment UNMAPPED=SampleID.step7-ConsensusReads-Unmapped.bam ALIGNED=SampleID.step9-Mapped-Consensus.bam R=$Ref O=SampleID.step10-Mapped-Consensus-UMI.bam

12- Step 11: Create index
    samtools index -b SampleID.step10-Mapped-Consensus-UMI.bam SampleID.step10-Mapped-Consensus-UMI.bai

13- Step 12: Call variant caller
    VarDict -G hg19.fa -N SampleID -b SampleID.step10-Mapped-Consensus-UMI.bam -f Allele_Freq -q 11 -r 2 -k 1 -I 50 -c 1 -S 2 -E 3 -g 4 PanelName.ampInsert.bed | teststrandbias.R | var2vcf_valid.pl -N $tumor_group -E -f $af > SampleID_AF_consensus.vardict.vcf



Usage
-----
```
#!text

usage: clonal_hematopoiesis.py [-h] [--help] [-c CONFIG [CONFIG ...]]
                               [-s STEPS] [-o OUTPUT_DIR]
                               [-j {pbs,batch,daemon,slurm}] [-f] [--no-json]
                               [--report] [--clean]
                               [-l {debug,info,warning,error,critical}]
                               [--sanity-check]
                               [--container {docker, singularity} {<CONTAINER PATH>, <CONTAINER NAME>}]
                               [-t {mugqic,mpileup,light}] [-r READSETS] [-v]

Version: 3.1.6-beta

For more documentation, visit our website: https://bitbucket.org/mugqic/mugqic_pipelines/

optional arguments:
  -h                    show this help message and exit
  --help                show detailed description of pipeline and steps
  -c CONFIG [CONFIG ...], --config CONFIG [CONFIG ...]
                        config INI-style list of files; config parameters are
                        overwritten based on files order
  -s STEPS, --steps STEPS
                        step range e.g. '1-5', '3,6,7', '2,4-8'
  -o OUTPUT_DIR, --output-dir OUTPUT_DIR
                        output directory (default: current)
  -j {pbs,batch,daemon,slurm}, --job-scheduler {pbs,batch,daemon,slurm}
                        job scheduler type (default: slurm)
  -f, --force           force creation of jobs even if up to date (default:
                        false)
  --no-json             do not create JSON file per analysed sample to track
                        the analysis status (default: false i.e. JSON file
                        will be created)
  --report              create 'pandoc' command to merge all job markdown
                        report files in the given step range into HTML, if
                        they exist; if --report is set, --job-scheduler,
                        --force, --clean options and job up-to-date status are
                        ignored (default: false)
  --clean               create 'rm' commands for all job removable files in
                        the given step range, if they exist; if --clean is
                        set, --job-scheduler, --force options and job up-to-
                        date status are ignored (default: false)
  -l {debug,info,warning,error,critical}, --log {debug,info,warning,error,critical}
                        log level (default: info)
  --sanity-check        run the pipeline in `sanity check mode` to verify that
                        all the input files needed for the pipeline to run are
                        available on the system (default: false)
  --container {docker, singularity} {<CONTAINER PATH>, <CONTAINER NAME>}
                        run pipeline inside a container providing a container
                        image path or accessible docker/singularity hub path
  -t {mugqic,mpileup,light}, --type {mugqic,mpileup,light}
                        DNAseq analysis type
  -r READSETS, --readsets READSETS
                        readset file
  -v, --version         show the version information and exit

Steps:
```
![clonal_hematopoiesis workflow diagram](https://bitbucket.org/mugqic/genpipes/raw/master/resources/workflows/GenPipes_clonal_hematopoiesis.resized.png)
[download full-size diagram](https://bitbucket.org/mugqic/genpipes/raw/master/resources/workflows/GenPipes_clonal_hematopoiesis.png)
```
------
1- picard_sam_to_fastq
2- fastqc
3- cutadapt
4- umi_extract
5- umi_deduplicate
6- sambamba_merge_sam_files
7- cram_output
8- metrics
9- picard_calculate_hs_metrics
10- vardict
11- snp_effect_vardict
12- run_multiqc

```
picard_sam_to_fastq
-------------------
Convert SAM/BAM files from the input readset file into FASTQ format
if FASTQ files are not already specified in the readset file. Do nothing otherwise.

fastqc
------
[FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/) provide quality 
control checks on raw sequence data coming from high throughput sequencing pipelines.
It provides a modular set of analyses including: per base sequence quality, per sequence 
quality scores, per base sequence content, per base GC content, per sequence GC content, 
per base N contenf, Sequence length distribitions, sequence duplication, onverrepresented
sequences and kmer content metrics.

cutadapt
--------
Raw reads quality trimming and removing adapters is performed using [Cutadapt](https://cutadapt.readthedocs.io/en/stable/index.html).
'Adapter1' and 'Adapter2' columns from the readset file ar given to Cutadapt. For PAIRED_END readsets, both adapters are used.
For SINGLE_END readsets, only Adapter1 is used and left unchanged.
To trim the front of the read use adapter_5p_fwd and adapter_5p_rev (for PE only) in cutadapt section of ini file.

This step takes as input files:

1. FASTQ files from the readset file if available
2. Else, FASTQ output files from previous picard_sam_to_fastq conversion of BAM files

umi_extract
-----------
1. Use fgbio ExtractUmisFromBam to Extract the UMI information and store as a tag

umi_deduplicate
---------------
1. Use fgbio GroupReadsByUmi to group reads that appeared to have come from the same original
molecule (this groups reads by template as well as by UMI sequence.)
2. Used fgbio CallDuplexConsensusReads to create duplex consensus reads by combining the 
evidence from a top and bottom single stranded consensus read

sambamba_merge_sam_files
------------------------
BAM readset files are merged into one file per sample. Merge is done using [Picard](http://broadinstitute.github.io/picard/).

This step takes as input files:

1. Aligned and sorted BAM output files from previous bwa_mem_picard_sort_sam step if available
2. Else, BAM files from the readset file

cram_output
-----------
Generate long term storage version of the final alignment files in CRAM format
Using this function will include the orginal final bam file into the  removable file list 

metrics
-------
Compute metrics and generate coverage tracks per sample. Multiple metrics are computed at this stage:
Number of raw reads, Number of filtered reads, Number of aligned reads, 
Median, mean and standard deviation of insert sizes of reads after alignment, percentage of bases
covered at X reads (%_bases_above_50 means the % of exons bases which have at least 50 reads)
whole genome or targeted percentage of bases covered at X reads (%_bases_above_50 means the % of exons
bases which have at least 50 reads). 

picard_calculate_hs_metrics
---------------------------
Compute on target percent of hybridisation based capture.

vardict
-------
Variant calls in single sample mode. [VarDict](https://github.com/AstraZeneca-NGS/VarDict) 
implements amplicon bias aware variant calling from targeted sequencing experiments, rescue 
of long indels by realigning bwa soft clipped reads and better scalability than many Java 
based variant callers effect annotation. 

snp_effect_vardict
------------------
Variant effect annotation. The .vcf files are annotated for variant effects using the SnpEff software.
SnpEff annotates and predicts the effects of variants on genes (such as amino acid changes).

run_multiqc
-----------

