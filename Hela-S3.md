# Pipeline for Hela-S3

## samples in this parts
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


The duplicates were removed using picard MarkDuplicates. 

Run following snakefile:
```python

keep=list(map(str, range(1,23)))+['X']
sm_list=[
'Hela_mock_rep1',
'Hela_mock_rep2',
'Hela_PMX_rep1',
'Hela_PMX_rep2',
'Hela_mock_rep1_inp',
'Hela_mock_rep1_inp',
'Hela_PMX_rep1_inp',
'Hela_PMX_rep1_inp'
]


rule project:
		input:
			expand("{sm}.merged.sort.bed.gz", sm=sm_list)


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
		"{sm}.merged.bed"
	shell:
		"cat {input} > {output}"

# sort bed
rule sortBed:
	input:
		"{sm}.merged.bed"
	output:
		"{sm}.merged.sort.bed.gz"
	shell:
		"bedtools sort -i {input} | bgzip -c >{output} && tabix {output}"

```

To view nucleotide frequencies around the predicted uracil sites:
```bash
python fetch_SeqContext.py sample.merged.sort.bed sample.seq.txt hg19.fa
```

To see how dU is distributed in different replication timing, [RT profile](https://www2.replicationdomain.com/) ( or see [10kb RT bedgraph](https://github.com/Jyyin333/Ucaps-seq/blob/main/files/RT_HeLaS3_hg19.10kb.bedgraph) )in Hela-S3 was used for investigating the relationship between dU and RT:
```python
# At first, the records with the thymine (T) nucleotide at uracil site were selected
sm_list=[
'Hela_mock_rep1',
'Hela_mock_rep2',
'Hela_PMX_rep1',
'Hela_PMX_rep2'
]

rule all:
	input:
		expand()

rule select_T:
	input:
		"{sm}.merged.sort.bed"
	output:
		bgz="{sm}.merged.sort.T.bed.gz",
		tbi="{sm}.merged.sort.T.bed.gz.tbi"
	shell:
		"python select_T.py {input} |bedtools sort -i - | bgzip -c >{output.bgz} && tabix {output.bgz}"


# calculate reads info in each RT state
rule getInfoinRegion:
	input:
		treat=rules.select_T.output.bgz,
		inp="{sm}_inp.merged.sort.bed.gz"
	output:
		"{sm}-RT_HeLaS3_hg19.10kb.tsv"
	shell:
		"python getInfoinRegion.py -i {input.treat} {input.inp} -b RT_HeLaS3_hg19.10kb.bedgraph -g hg19.2bit -o {output}"


```
the results sample-RT_HeLaS3_hg19.10kb.tsv looks like this:
```
chr1	26499	36499	-0.12170513999999999	0	0	1	0	2565	2617
chr1	36499	46499	0.07242317000000001	0	0	0	0	2695	3385
chr1	46499	56499	0.2092264	1	2	3	4	3224	3038
chr1	56499	66499	0.2307322	2	1	3	2	3561	3105
chr1	66499	76499	0.3413856	1	0	0	4	3205	3193
```
The last six columns are Treat-plus-readCounts, Treat-minus-readCounts, Inp-plus-readCounts, Inp-minus-readCounts, base-A-counts, base-T-counts.

Finally, the boxplots are drawn and the differences among RT-states are statistically tested:
```bash
for sm in `ls *RT_HeLaS3_hg19.10kb.tsv`
do
	Rscript RT_test.r --input $sm --output ${sm%.*}.wtest_res.txt --outfig ${sm%.*}.pdf
done

```
