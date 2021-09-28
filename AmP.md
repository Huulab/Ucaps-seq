# Pipeline for AmP data

## samples in this part
| sample title | pipeline |
| :---: | :---: |
| AmP rnf rep1 | AmP |
| AmP rnf rep2 | AmP |
| AmP site3 rep1 | AmP |
| AmP site3 rep2 | AmP |
| AmP OT1 rep1 | AmP |
| AmP WT rnf rep1 | AmP |
| AmP WT site3 rep1 | AmP |
| AmP WT OT1 rep1 | AmP |



Run snakefile:
```python
sm_list=[
'AmP_rnf_rep1',
'AmP_rnf_rep2',
'AmP_site3_rep1',
'AmP_site3_rep2',
'AmP_OT1_rep1',
'AmP_WT_rnf_rep1',
'AmP_WT_site3_rep1',
'AmP_WT_OT1_rep1'
]

rule all:
	input:
		expand("{sm}.done",sm=smS)


rule retag:
	input:
		r1="{sm}_combined_R1.fastq.gz",
		r2="{sm}_combined_R2.fastq.gz"
	output:
		"{sm}_tagged_R1.fastq.gz",
		"{sm}_tagged_R2.fastq.gz"
	params:
		prefix="{sm}"
	shell:
		"python reTagFastq.py -i {input.r1} {input.r2} -p {params.prefix}"

rule bwa_map:
	input:
		"{sm}_tagged_R1.fastq.gz",
		"{sm}_tagged_R2.fastq.gz"
	output:
		temp("{sm}.bam")
	shell:
		"bwa mem {bwa_index_path} {input} | samtools view -Sb - > {output}"

rule sort_bam:
	input:
		rules.bwa_map.output
	output:
		bam="{sm}_sorted.bam",
		bai="{sm}_sorted.bai"
	shell:
		"picard SortSam"
		" -I {input}"
		" -O {output.bam}"
		" --SORT_ORDER coordinate"
		" --CREATE_INDEX true"

rule filterbam:
	input:
		bam="{sm}_sorted_dedup.bam",
		bai="{sm}_sorted_dedup.bai",
	output:
		b="{sm}.clean.bam",
		i="{sm}.clean.bam.bai",
	shell:
		"python filterBam.py -b {input.bam} -o {output.b}"

rule DONE:
	input:
		rules.sort_bam.output
	output:
		touch("{sm}.done")

```

Note, the duplicates in targeted amplicon sequencing data were removed using barcode. If two reads have same genomic coordinates and same paired barcodes, they were marked as duplicates and would be discarded.
```bash
# remove duplicates
# region referes to targeted region
# primer refers to primer seq used for PCR
python dedupBarcode.py --bam AmP_rnf_rep1.clean.bam --region chr1:185056690-185057100 --primer TCAGGCTGTGCAGACAAACGG TACAGCAAGGAGGACTTGCCC --out rnf_targets.bam


# get mutant infomation
python getMutInfo.py --bam rnf_targets.bam --out rnf_targets.mut.info

# merge mutant infomation
python mergeMut.py -i rnf_targets.mut.info -o rnf_targets.mergeMut.info --bam rnf_targets.bam
```

The final results:
```
# chrom    0_base_position       ref_base       query_base       counts    strand       strand_specific_counts  all_counts
chr1    185056689       C       T       1545    +       117501  117501
chr1    185057104       C       T       525     -       6795    6795
chr1    185056799       C       T       141     +       391094  391094
chr1    185056808       C       T       89      +       391046  391046
chr1    185056726       C       T       81      +       391215  391215
chr1    185056791       C       T       78      +       391097  391097
chr1    185056760       C       T       69      +       391126  391126
chr1    185056757       C       T       59      +       391126  391126

...
```

This file is used to calculate mutation rates.
