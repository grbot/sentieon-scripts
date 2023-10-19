# sentieon-scripts

## Workflows

### dnascope-align-call

1) Set sample sheet. See `NA12878.samplesheet.small.tsv` 
2) Set workflow in `nextflow.config`. `workflow = 'dnascope-align-call'`
3) Run: `nextflow run main.nf -c nextflow.config -w /cbio/projects/020/gerrit/sentieon-scripts-align-work -profile ilifu -resume`

