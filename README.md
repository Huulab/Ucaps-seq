# Ucaps-seq
This repository includes source codes used in ***Sequencing uracil in DNA at single-nucleotide resolution***.

## Data availability
All raw sequencing data are available at NCBI Sequence Read Archive with BioProject ID [PRJNA728500](https://www.ncbi.nlm.nih.gov/bioproject/?term=PRJNA728500)

## Samples

| sample title | pipeline |
| :---: | :---: |
| Hela UNGKO Mock rep1 | Hela-S3 |
| Hela UNGKO Mock rep2 | Hela-S3 |
| Hela UNGKO PMX rep1 | Hela-S3 |
| Hela UNGKO PMX rep2 | Hela-S3 |
| Hela UNGKO Mock inp rep1 | Hela-S3 |
| Hela UNGKO Mock inp rep2 | Hela-S3 |
| Hela UNGKO PMX inp rep1 | Hela-S3 |
| Hela UNGKO PMX inp rep2 | Hela-S3 |
| HEK293T WT rep1 | HEK293T |
| HEK293T WT rep2 | HEK293T |
| HEK293T NOD rep1 | HEK293T |
| HEK293T NOD rep2 | HEK293T |
| HEK293T 10X rep1 | HEK293T |
| HEK293T 10X rep2 | HEK293T |
| HEK293T 50X rep1 | HEK293T |
| HEK293T 50X rep2 | HEK293T |
| HEK293T 100X rep1 | HEK293T |
| HEK293T 100X rep2 | HEK293T |
| MM Bcell AK rep1 | CH12F3 |
| MM Bcell AK rep2 | CH12F3 |
| MM Bcell UK rep1 | CH12F3 |
| MM Bcell UK rep2 | CH12F3 |
| Oligo dU | synthetic DNA |
| Oligo input | synthetic DNA |
| AmP rnf rep1 | AmP |
| AmP rnf rep2 | AmP |
| AmP site3 rep1 | AmP |
| AmP site3 rep2 | AmP |
| AmP OT1 rep1 | AmP |
| AmP WT rnf rep1 | AmP |
| AmP WT site3 rep1 | AmP |
| AmP WT OT1 rep1 | AmP |


## pre-analysis
First of all, FastQC was used for quality control and BWA mem was used for aligning reads to appropriate refernece genome.
```bash
# quickly check
fastqc --noextract --format=fastq *.fastq.gz

# mapping and sort
bwa mem $bwa_index_path sample.R1.fastq.gz sample.R2.fastq.gz | samtools view -Sb - > sample.bam
picard SortSam -I sample.bam -O sample_sorted.bam --SORT_ORDER true --CREATE_INDEX true 
```
When this step is complete, the downstream analysis can be performed based on pipeline to which the sample belongs.


## pipeline
- [Ucaps-seq for synthetic DNA samples](https://github.com/Jyyin333/Ucaps-seq/blob/main/sDNA.md)
- [Ucaps-seq for Hela-S3 samples](https://github.com/Jyyin333/Ucaps-seq/blob/main/Hela-S3.md)
- [Ucaps-seq for CH12F3 samples](https://github.com/Jyyin333/Ucaps-seq/blob/main/CH12F3.md)
- [Ucaps-seq for HEK293T samples](https://github.com/Jyyin333/Ucaps-seq/blob/main/HEK293T.md)
