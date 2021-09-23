# Pipeline for synthetic DNA

## samples in this parts
| sample title | pipeline |
| :---: | :---: |
| Oligo dU | synthetic DNA |
| Oligo input | synthetic DNA |

[synthetic DNA reference genome](https://github.com/Jyyin333/Ucaps-seq/blob/main/files/sDNA.fa) was used for creating BWA index.
```bash
bwa index sDNA.fa
```

After mapping, using Python script:
```bash
# for dU sample
python Treat.stat.py --bam dU.bam --out dU.res --json dU.json

# for input
python Inp.stat.py --bam inp.bam --out inp.res --json inp.json
```

At a glance of the results:
```
# dU/Inp res
dU-forward:1861344	dU-reverse:8305
dT-forward:10787	dT-reverse:4610

# dU json
{
    "11S": {
        "11S99M40S": 42551
    },
    "12S": {
        "12S99M39S": 68839
    },
    "13S": {
        "13S99M38S": 86302
    },
    "14S": {
        "14S99M37S": 1587337
    },
    "15S": {
        "15S99M36S": 42701
    },
    "16S": {
        "16S99M35S": 954
    },
    "17S": {
        "17S99M34S": 118
    }
}
```

In this section, the results showed that compared with input, the dU-containing strands (UF) were highly enriched in Ucaps-seq, 
and the polymerase mainly stopped right before the uracil site
