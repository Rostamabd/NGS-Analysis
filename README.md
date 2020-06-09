# Next Generation Sequencing Work-flow of Bos taurus genome

**Author**: *Rostam Abdollahi-Arpanahi*

**Date**: May 7th, 2020

-------------------

Make sure all the following tools are already installed on your machine.

```
module load fastqc
module load bwa
module load samtools
module load picard
module load gatk
```

## 1. Data Preparation

Basically the NGS company or databases provide the data in fasta or SRA format. You must be able to download the NGS to your own computer. The most common command for downloading the data in unix is:

```
wget  "URL address"
```



Combine multiple *.gz files without un-compressing (Here we had to lanes of reads per each sample)

- Forward strand

```
zcat read1_1.fastq.gz  read1_2.fastq.gz | gzip -c > forward.fastq.gz
zcat forward.fastq.gz | md5sum
```

- Reverse strand

```
zcat read2_1.fastq.gz  read2_2.fastq.gz | gzip -c > reverse.fastq.gz
zcat reverse.fastq.gz | md5sum
```

##  2. Quality Control       ##
- check the quality control of NGS data

```
fastqc -t 4 forward.fastq.gz reverse.fastq.gz
```

- If the data quality is not good enough, use a program such as trimmometric or cutadapt or trim_galore (we skipped this step in our analysis because of great data quality)

```
java -jar trimmomatic-0.36.jar PE LLUMINACLIP:TruSeq3-PE.fa:2:30:10 SLIDINGWINDOW:6:20 MINLEN:140 CROP:140 -o output_files file1.fatq.gz file2.fastq.gz file3.fastq.gz file4.fastq.gz
```

##  3. Alignment/Mapping     ##


- After quality control, the next step is to align the reads to a reference sequence

### 3.1 Download the latest Bovine reference genome
- It is available for download at NCBI through the following link: [ARS-UCD1.2 refernce genome](https://www.animalgenome.org/repository/cattle/ARS-UCD1.2_HISAT2Index/ARS-UCD.1.2.fa.gz)

```
gunzip ARS-UCD.1.2.fa.gz
mv ARS-UCD.1.2.fa bos.fa
```

### 3.2 - Index the reference file
```
bwa index -a bwtsw bos.fa
```

### 3.3 - Create a dictionary and index for reference genome
```
picard CreateSequenceDictionary -R bos.fa
samtools faidx bos.fa
```

### 3.4 - Generate alignments in the suffix array coordinate
```
bwa aln -t 4 bos.fa   forward.fastq.gz > forward_read.sai
bwa aln -t 4 bos.fa   reverse.fastq.gz > reverse_read.sai
```

### 3.5 -  Map to the reference genome with BWA

```
bwa mem bos.fa forward_read.sai reverse_read.sai forward.fastq.gz  reverse.fastq.gz > aln_read.sam
```

##  4. Alignment Post-Processing       ##
###  4.1 - Convert sam file to bam file

```
samtools view -bS  -t bos.fa aln_read.sam  -o aln_read.bam
```

### 4.2 - sort and mapping coverage and mapping depth
```
samtools sort aln_read.bam -o sorted_reads.bam
```

### 4.3 - Check the number of reads and coverage
```
samtools flagstat sorted_reads.bam > coverage_mapping
```

### 4.4 -  Marking duplicates with picard

```
picard  MarkDuplicates  INPUT=sorted_reads.bam  OUTPUT=dedup_reads.bam METRICS_FILE=output.Metrics VALIDATION_STRINGENCY=LENIENT  MAX_FILE_HANDLES_FOR_READ_ENDS_MAP=50
```

### 4.5 - Fixed the bam file (just in case there is any error)

```
picard FixMateInformation  I=dedup_reads.bam  O=nmsrt_fixmate_sorted_reads.bam SO=queryname VALIDATION_STRINGENCY=LENIENT TMP_DIR=./tmp/temp
```

### 4.6 - Sort the bam file again
```
picard SortSam I=dedup_reads.bam O=picard_sorted.bam SORT_ORDER=coordinate TMP_DIR=./tmp/temp  VALIDATION_STRINGENCY=LENIENT
```

### 4.7 - Add the read groups (just in case it is missing)
```
picard AddOrReplaceReadGroups  I=picard_sorted.bam  O=sorted_reads_RG_Edited.bam RGID=4 RGLB=lib1 RGPL=illumina RGPU=unit1 RGSM=chr8 VALIDATION_STRINGENCY=LENIENT
mv sorted_reads_RG_Edited.bam picard_sorted_RG.bam
```

### 4.8 - Index using picard
```
picard BuildBamIndex   I=picard_sorted_RG.bam  VALIDATION_STRINGENCY=LENIENT
```

### 4.9 - Reorder
```
java -jar picard.jar ReorderSam I=picard_sorted_RG.bam O=reorder.bam SD=bos.fa CREATE_INDEX=TRUE TMP_DIR=./tmp/temp

mv  reorder.bam picard_sorted_RG.bam
```


### 4.9 - Download annotation file from database (Ensembl or [1000 Bull Genomes Project](http://www.1000bullgenomes.com/))
```
wget http://www.1000bullgenomes.com/doco/ARS1.2PlusY_BQSR_v3.vcf.gz
mv ../sorted_bos.vcf.idx ../bos.vcf.idx
```


### 4.10 - Realignment with gatk. It is recommendedto use the GATK version 3.8 for this purpose

```
java -jar GenomeAnalysisTK.jar -R bos.fa -T RealignerTargetCreator -o indels_Realigner.intervals   -I picard_sorted_RG.bam -known bos.vcf

java -jar GenomeAnalysisTK.jar -R bos.fa -T IndelRealigner -targetIntervals indels_Realigner.intervals -known bos.vcf  -I picard_sorted_RG.bam -o  sorted_reads_indel.bam
```

#### 4.11 - Recalibration process

```
java -Xmx16g -jar  GenomeAnalysisTK.jar -R bos.fa -T BaseRecalibrator  -I sorted_reads_indel.bam -knownSites  bos.vcf  -o sorted_reads_indel_grp

java -Xmx16g -jar GenomeAnalysisTK.jar -R bos.fa -T PrintReads -BQSR sorted_reads_indel_grp -I  picard_sorted_RG.bam -o sorted_reads_indel_recal.bam 
```



##  5 - Variant calling       	##

### 5.1 - Variant Calling with gvcf format (we assume we have two different *.bam files, you may have several bam files)

```
java -Xmx16g -jar  GenomeAnalysisTK.jar -R bos.fa -T HaplotypeCaller  -I sorted_reads_indel_recal.bam -I sorted_reads_ctrl_indel_recal.bam --dbsnp bos.vcf -stand_call_conf 20  -o tothap.raw.snps.indels.g.vcf
```

### 5.2 - SNP (SNV) calling

```
java -jar GenomeAnalysisTK.jar -R bos.fa  -T SelectVariants -selectType SNP -V tothap.raw.snps.indels.g.vcf -o SNP.vcf
```

### 5.3 - InDel calling
```
java -jar GenomeAnalysisTK.jar -R bos.fa  -T SelectVariants -selectType INDEL   -V tothap.raw.snps.indels.g.vcf -o Indel.vcf
```

### 5.4 - Quality Control of SNPs (filter SNVs)

```
java -jar GenomeAnalysisTK.jar  -T VariantFiltration -R bos.fa --variant SNP.vcf   --filterName "snpsfilter" --filterExpression "QD<2.0||MQ<40.0||FS>60.0" --out snps.tagged.vcf 
```

### 5.5 - Quality Control of InDels (filter InDels)
```
java -jar GenomeAnalysisTK.jar  -T VariantFiltration -R bos.fa --variant Indel.vcf   --filterName "my_indel_filter" --filterExpression "QD<2.0||FS>200.0||ReadPosRankSum < -20.0" --out indels.tagged.vcf
```

### 5.6 select filtered variants

```
java -jar GenomeAnalysisTK.jar  -T SelectVariants  -R bos.fa  --variant snps.tagged.vcf -select 'vc.isNotFiltered()' -o snps.filtered.vcf

java -jar GenomeAnalysisTK.jar  -T SelectVariants  -R bos.fa  --variant indels.tagged.vcf  -select 'vc.isNotFiltered()' -o indels.filtered.vcf
```

## Contact Information

Please send your comments and suggestions to rostam7474 at gmail dot com
