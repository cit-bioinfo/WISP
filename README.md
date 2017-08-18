# SMAP
SMAP: sequencing patient-derived xenografts



SMAP is a pipeline designed to process xenografts sequencing data. 

It takes FASTQ as input and outputs specie-specific BAM or gene counts.



## SMAP steps

These are the few steps to run a complete SMAP pipeline, from genome file download to obtaining specie-specific gene counts and fusion transcripts.

1. [Download and combine reference genomes and transcriptomes][#downloadCombine]. A shell script is given to simplify the genome combination process.
2. [Alignment to combined reference genome/transcriptome][## Align reads to the combined reference genome and transcriptome using STAR]. This should follow the usual process of read alignement.
3. ​



##  Prepare combined reference genomes and transcriptome

<a name="downloadCombine"></a>

These steps, shown for human and mouse xenografts, show how to prepare reference genome and transcriptome for SMAP. 



### Download reference genomes and transcriptomes

It may take several minutes to download the entire genomes of both species.

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

### Combine genomes and transcriptomes

```shell
sh SMAP_prepareReference.sh \
-t Homo_sapiens.GRCh37.75.dna_sm.primary_assembly.fa \
-u Homo_sapiens.GRCh37.75.gtf \
-m Mus_musculus.GRCm38.75.dna_sm.primary_assembly.fa \
-s Mus_musculus.GRCm38.75.gtf  \
-o combined
```

The combined genome and transcriptome will be available under the **combined** directory, as announced by the output of `SMAP_prepareReference.sh`.



## Align reads to the combined reference genome and transcriptome using STAR

The  `SMAP_prepareReference.sh` script gives the next command to run which aims at preparing files for RNAseq reads alignement using STAR. This only needs to be done once for all samples.

```shell
mkdir combined/StarIndexed_combined

STAR --runMode genomeGenerate --genomeDir combined/StarIndexed_combined    --genomeFastaFiles combined/combined.fasta  --sjdbGTFfile combined/combined.gtf  --sjdbOverhang 100   # set to the length of base pair sequencing 
```



Each of the samples can then be mapped using the following STAR command in which `combined/StarIndexed_combined` is the STAR indexed combined genome/transcriptome generated in the precedent command, `FASTQ1.fq` and `FASTQ2.fq` are the paired raw FASTQ files of the sample and `combined/combined.gtf` is the newly  reference transcriptome (gtf) created by the `SMAP_prepareReference.sh` script.

```shell
STAR --genomeDir combined/StarIndexed_combined --readFilesIn FASTQ1.fq FASTQ2.fq --sjdbGTFfile combined/combined.gtf --outFilterType BySJout --outFilterMultimapNmax 100   --outSAMtype BAM SortedByCoordinate  

# Other parameters for STAR may include (not specific to SMAP) :
# --runThreadN 16 
# --alignSJoverhangMin 8 --alignSJDBoverhangMin 1
# --outFilterMismatchNmax 999 --outFilterMismatchNoverLmax 0.04 --alignIntronMin 20 
# --alignIntronMax 1000000 --alignMatesGapMax 0 --outSAMtype BAM SortedByCoordinate
# --chimSegmentMin 5 --chimJunctionOverhangMin 5  --chimScoreMin 0
```



The output is a classical BAM file, that should be named `Aligned.sortedByCoord.out.bam`.



## Gene expression using SMAP and FeatureCount




## Fusion


refCDNA=/datacit/00_DATABANKS/ensembl_75_humanMouseXenome/hs75.hg19_mmu75.GRCm38_chrename_CDNA.fasta
```sh

wget -O -   | gunzip -c ftp://ftp.ensembl.org/pub/release-75/fasta/homo_sapiens/cdna/Homo_sapiens.GRCh37.75.cdna.all.fa.gz >Homo_sapiens.GRCh37.75.cdna.all.fa

wget -O -   | gunzip -c  ftp://ftp.ensembl.org/pub/release-75/fasta/mus_musculus/cdna/Mus_musculus.GRCm38.75.cdna.all.fa.gz >Mus_musculus.GRCm38.75.dna_sm.primary_assembly.fa

cat Homo_sapiens.GRCh37.75.cdna.all.fa Mus_musculus.GRCm38.75.dna_sm.primary_assembly.fa >combined/combinedCdna.fa

```

```sh


STAR-Fusion -J Chimeric.out.junction -G  combined/combined.gtf -C combined/combinedCdna.fa

```
