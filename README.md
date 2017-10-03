# SMAP

SMAP is a pipeline for the process of xenografts sequencing data.

It takes FASTQ as input and outputs specie-specific BAM or gene counts.


## SMAP steps


These are the few steps to run a complete SMAP pipeline, from genome file download to obtaining specie-specific gene counts and fusion transcripts.

0. SMAP requires a simple [installation](#install) step
1. [Download and combine reference genomes and transcriptomes](#downloadCombine). A shell script is given to simplify the genome combination process
2. [Alignment to combined reference genome/transcriptome](#align). This should follow the usual process of read alignement
3. [Gene expression quantification](#genexp). The SMAP R package provides all necessary steps to quantify the expression of genes (or exons) of the host and the graft sequences seperatly
4. [Fusion detection](#fusion). The SMAP R package provides a function to efficiently filter false positive gene fusions
5. [BAM file seperation](#bamsplit), only required for other types of analysis such as variant calling.


## 0 Install
<a name="install"></a>


The SMAP R package can be installed by downloading/cloning [the repository](https://github.com/RemyNicolle/SMAP) and using `R CMD INSTALL` or directly from R using the devtools package: 
```R
# If not already installed, install the devtools package
install.packages("devtools")

library(devtools)
install_github("RemyNicolle/SMAP")
```

The SMAP pipeline requires two additional scripts. These are included in the downloaded/cloned repository. They can also be downloaded seperately:
- A shell script to prepare the genomes and transcriptomes of reference: [SMAP_prepareReference.sh](https://raw.githubusercontent.com/RemyNicolle/SMAP/master/SMAP_prepareReference.sh)
- A python script to seperate the BAM file of aligned reads based on their predicted specie of origin: [smap_splitBySpecie_standard.py](https://raw.githubusercontent.com/RemyNicolle/SMAP/master/smap_splitBySpecie_standard.py)


To download these scripts seperately:
```shell
wget https://raw.githubusercontent.com/RemyNicolle/SMAP/master/SMAP_prepareReference.sh
wget https://raw.githubusercontent.com/RemyNicolle/SMAP/master/smap_splitBySpecie_standard.py
```

__Requirements__ The BAM file seperation `smap_splitBySpecie_standard.py` script (which is not required for gene quantification) requires the `pysam` python module. 



## 1 Prepare combined reference genomes and transcriptome

<a name="downloadCombine"></a>

These steps, shown for human and mouse xenografts, detail how to prepare reference genome and transcriptome for SMAP. 



### 1.1 Download reference genomes and transcriptomes

It may take several minutes to download the entire genome of both species.

Human :

```shell
# Genome Fasta :
wget -O - ftp://ftp.ensembl.org/pub/release-75/fasta/homo_sapiens/dna/Homo_sapiens.GRCh37.75.dna_sm.primary_assembly.fa.gz  | gunzip -c > Homo_sapiens.GRCh37.75.dna_sm.primary_assembly.fa  

# Transcriptome GTF :
wget -O - ftp://ftp.ensembl.org/pub/release-75/gtf/homo_sapiens/Homo_sapiens.GRCh37.75.gtf.gz | gunzip -c >Homo_sapiens.GRCh37.75.gtf

```

Mouse :

```shell
# Genome Fasta :  
wget -O -  ftp://ftp.ensembl.org/pub/release-75/fasta/mus_musculus/dna/Mus_musculus.GRCm38.75.dna_sm.primary_assembly.fa.gz | gunzip -c   >Mus_musculus.GRCm38.75.dna_sm.primary_assembly.fa
# Transcriptome GTF :
wget -O -   ftp://ftp.ensembl.org/pub/release-75/gtf/mus_musculus/Mus_musculus.GRCm38.75.gtf.gz | gunzip -c > Mus_musculus.GRCm38.75.gtf   
```

### 1.2 Combine genomes and transcriptomes


Use the provided `SMAP_prepareReference.sh` shell script to prepare the reference genome and transcriptomes.


```shell
sh SMAP_prepareReference.sh \
-t Homo_sapiens.GRCh37.75.dna_sm.primary_assembly.fa \
-u Homo_sapiens.GRCh37.75.gtf \
-m Mus_musculus.GRCm38.75.dna_sm.primary_assembly.fa \
-s Mus_musculus.GRCm38.75.gtf  \
-o combined
```

The combined genome and transcriptome will be available under the **combined** directory (can be modified using the `-o` argument), as announced by the output of `SMAP_prepareReference.sh`.



## 2 Align reads to the combined reference genome and transcriptome using STAR
<a name="align"></a>

The SMAP pipeline is aligner-agnostic. However, STAR is recommended. This section details how to run STAR. Any personnal or institutional STAR pipeline can be used here with the SMAP combined genome/transcriptome.


STAR first requires the combined genome/transcriptome to be indexed. This only needs to be done once for a given pair of genomes.

```shell
mkdir combined/StarIndexed_combined

STAR --runMode genomeGenerate --genomeDir combined/StarIndexed_combined    --genomeFastaFiles combined/combined.fasta  --sjdbGTFfile combined/combined.gtf  --sjdbOverhang 100   # set to the length of base pair sequencing 
```



Each of the samples can then be mapped using the following STAR command in which `combined/StarIndexed_combined` is the STAR indexed SMAP combined genome/transcriptome generated in the precedent command, `FASTQ1.fq` and `FASTQ2.fq` are the paired raw FASTQ files of the sample and `combined/combined.gtf` is the newly  reference transcriptome (gtf) created by the `SMAP_prepareReference.sh` script.

```shell
STAR --genomeDir combined/StarIndexed_combined --readFilesIn FASTQ1.fq FASTQ2.fq --sjdbGTFfile combined/combined.gtf --outFilterType BySJout --outFilterMultimapNmax 100   --outSAMtype BAM SortedByCoordinate  

# Other parameters for STAR may include (not specific to SMAP) :
# --runThreadN 16 
# --alignSJoverhangMin 8 --alignSJDBoverhangMin 1
# --outFilterMismatchNmax 999 --outFilterMismatchNoverLmax 0.04 --alignIntronMin 20 
# --alignIntronMax 1000000 --alignMatesGapMax 0 --outSAMtype BAM SortedByCoordinate
# --chimSegmentMin 5 --chimJunctionOverhangMin 5  --chimScoreMin 0
```



The output is a classical BAM file, that is usually be named `Aligned.sortedByCoord.out.bam`.



## 3 Gene expression using SMAP and FeatureCount
<a name="genexp"></a>
Use the `SMAPcount` function in the SMAP R package. This function gives two separate count data matrix for the tumor and host respectively from bam files using featurecount function at the gene or exon level. 

```R
 d=system.file("extdata","Example_STAR_OUTPUT", package = "SMAP")
    counts = SMAPcount(SMAPBAM=d, GTF =paste(d, "combined/combined.gtf", sep="/"))
    countsTumor= counts$FCcountsTumor
    countsHost= counts$FCcountsHost
```


## 4 Fusion
<a name="fusion"></a>

refCDNA=/datacit/00_DATABANKS/ensembl_75_humanMouseXenome/hs75.hg19_mmu75.GRCm38_chrename_CDNA.fasta
```sh

wget -O -   | gunzip -c ftp://ftp.ensembl.org/pub/release-75/fasta/homo_sapiens/cdna/Homo_sapiens.GRCh37.75.cdna.all.fa.gz >Homo_sapiens.GRCh37.75.cdna.all.fa

wget -O -   | gunzip -c  ftp://ftp.ensembl.org/pub/release-75/fasta/mus_musculus/cdna/Mus_musculus.GRCm38.75.cdna.all.fa.gz >Mus_musculus.GRCm38.75.dna_sm.primary_assembly.fa

cat Homo_sapiens.GRCh37.75.cdna.all.fa Mus_musculus.GRCm38.75.dna_sm.primary_assembly.fa >combined/combinedCdna.fa

```

```sh


STAR-Fusion -J Chimeric.out.junction -G  combined/combined.gtf -C combined/combinedCdna.fa

```

## 5 Seperate BAM files
<a name="bamsplit"></a>
