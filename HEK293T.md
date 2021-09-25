# Pipeline for HEK293T

## samples in this part
| sample title | pipeline |
| :---: | :---: |
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


The duplicates were removed using picard MarkDuplicates. 

Run following snakefile:
```python

keep=list(map(str, range(1,23)))+['X']
sm_list=[
'HEK293T_WT_rep1',
'HEK293T_WT_rep2',
'HEK293T_NOD_rep1',
'HEK293T_NOD_rep2',
'HEK293T_10X_rep1',
'HEK293T_10X_rep2',
'HEK293T_50X_rep1',
'HEK293T_50X_rep2',
'HEK293T_100X_rep1',
'HEK293T_100X_rep2'
]


rule project:
	input:
		expand("{sm}.merged.sort.bed.gz.tbi", sm=sm_list),
		expand("{sm}.merged.sort.bam", sm=sm_list)


# remove duplicates
rule dedup:
	input:
		bam="{sm}_sorted.bam",
		bai="{sm}_sorted.bai"
	output:
		bam="{sm}_sorted_dedup.bam",
		bai="{sm}_sorted_dedup.bai",
		metrics="{sm}_sorted_dedup.metrics"
	threads: 5
	shell:
		"picard MarkDuplicates \
			--INPUT {input.bam} \
			--OUTPUT {output.bam}  \
			--METRICS_FILE {output.metrics} \
			--REMOVE_DUPLICATES true \
			--CREATE_INDEX true \
			--VALIDATION_STRINGENCY SILENT \
		"


# futher filter & transfer format
rule fetch_dU_by_chrom:
	input:
		bam="{sm}_sorted_dedup.bam",
		bai="{sm}_sorted_dedup.bai"
	output:
		bed=temp("{sm}.chr{i}.raw.bed")
	params:
		chrom="chr{i}"
	shell:
		"python fetch_dU_by_chrom.py \
		--ref_fa {hg19.fa} \
		--bam {input.bam} \
		--out_file {output.bed} \
		--chrom {params.chrom} \
		"

# merge sub-bed files
rule merge_bed:
	input:
		["{sm}.chr"+ i + ".raw.bed" for i in keep]
	output:
		temp("{sm}.merged.bed")
	shell:
		"cat {input} > {output}"

# sort bed
rule bedSort:
	input:
		rules.merge_bed.output
	output:
		temp("{sm}.merged.sort.bed")
	shell:
		"bedSort {input} {output}"

# bed2bam
rule bed2bam:
	input:
		rules.bedSort.output
	output:
		"{sm}.merged.sort.bam"
	shell:
		"bedtools bedtobam -i {input} -g hg19.chrom.sizes > {output} && samtools index {output}"

# gzip file
rule bgzip:
	input:
		rules.bedSort.output
	output:
		"{sm}.merged.sort.bed.gz"
	shell:
		"bgzip -c {input} > {output}"

# create index
rule tabix:
	input:
		rules.bgzip.output
	output:
		"{sm}.merged.sort.bed.gz.tbi"
	shell:
		"tabix {input}"

```

To detect potential offtargets, we first pileup all reads and then filter sites abide by rules as shown in  FigureS5g, at last, DESeq2 was used for 
testing whether the difference was statistically significant.
```bash
# pileup
python plreads.py -t HEK293T_NOD_rep1.merged.sort.bed.gz HEK293T_NOD_rep2.merged.sort.bed.gz -c HEK293T_WT_rep1.merged.sort.bed.gz HEK293T_WT_rep2.merged.sort.bed.gz \
-g hg19.fa -o HEK293T_NOD.pl

# filter sites
python filtersites.py HEK293T_NOD.pl HEK293T_NOD.filter.sites

# call offtargets
Rscript call_targets.r --treat HEK293T_NOD_rep1.merged.sort.bam HEK293T_NOD_rep2.merged.sort.bam --control HEK293T_WT_rep1.merged.sort.bam HEK293T_WT_rep2.merged.sort.bam \
--window 30 --sites HEK293T_NOD.filter.sites --out HEK293T_NOD.candidate.sites

```
