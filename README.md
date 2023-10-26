# sentieon-scripts

## Workflows

### dnascope-align-call

1) Set sample sheet. See `NA12878.samplesheet.small.tsv` 
2) Set workflow in `nextflow.config`. `workflow = 'dnascope-align-call'`
3) Run: `nextflow run main.nf -c nextflow.config -w /cbio/projects/020/gerrit/sentieon-scripts-align-work -profile ilifu -resume`

### dnaseq-align-call

1) Set sample sheet. See `NA12878.samplesheet.small.tsv` 
2) Set workflow in `nextflow.config`. `workflow = 'dnaseq-align-call'`
3) Run: `nextflow run main.nf -c nextflow.config -w /cbio/projects/020/gerrit/sentieon-scripts-align-work -profile ilifu -resume`

### dnascope-genotype-gvcfs

1) Set `gvcf.list`
2) Set workflow in `nextflow.config`. `workflow = 'dnascope-genotype-gvcfs'`
3) Run: `nextflow run main.nf -c nextflow.config -w /cbio/projects/020/gerrit/sentieon-scripts-align-work -profile ilifu -resume`

### dnaseq-genotype-gvcfs

1) Set `gvcf.list`
2) Set workflow in `nextflow.config`. `workflow = 'dnaseq-genotype-gvcfs'`
3) Run: `nextflow run main.nf -c nextflow.config -w /cbio/projects/020/gerrit/sentieon-scripts-align-work -profile ilifu -resume`

