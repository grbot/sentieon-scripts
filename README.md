# sentieon-scripts

## Workflows

### dnascope-align-call (aligment from Fastq and calling from CRAMs)

1) Set sample sheet. See `NA12878.samplesheet.small.tsv` (contains Fastq reads)
2) Set workflow in `nextflow.config`. `workflow = 'dnascope-align-call'`
3) Run: `nextflow run main.nf -c nextflow.config -w /cbio/projects/020/gerrit/work -profile ilifu -resume`

### dnascope-call (calling from CRAM/BAMs)

1) Set sample sheet. See `NA12878.samplesheet.large.crams.tsv` (contains CRAM/BAM alignments)
2) Set workflow in `nextflow.config`. `workflow = 'dnascope-call'`
3) Run: `nextflow run main.nf -c nextflow.config -w /cbio/projects/020/gerrit/work -profile ilifu -resume`

### dnascope-genotype-gvcfs (joint genotyping)

1) Set `gvcf.list` (needs to be a concatenated list of gVCFS produced earlier)
2) Set workflow in `nextflow.config`. `workflow = 'dnascope-genotype-gvcfs'`
3) Run: `nextflow run main.nf -c nextflow.config -w /cbio/projects/020/gerrit/work -profile ilifu -resume`

### dnaseq-align-call (aligment from Fastq and calling from CRAMs)

1) Set sample sheet. See `NA12878.samplesheet.small.tsv` (contains Fastq reads)
2) Set workflow in `nextflow.config`. `workflow = 'dnaseq-align-call'`
3) Run: `nextflow run main.nf -c nextflow.config -w /cbio/projects/020/gerrit/work -profile ilifu -resume`

### dnaseq-genotype-gvcfs (joint genotyping)

1) Set `gvcf.list` (needs to be a concatenated list of gVCFS produced earlier)
2) Set workflow in `nextflow.config`. `workflow = 'dnaseq-genotype-gvcfs'`
3) Run: `nextflow run main.nf -c nextflow.config -w /cbio/projects/020/gerrit/work -profile ilifu -resume`

