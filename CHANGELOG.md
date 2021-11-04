## 1.0.2 - 2021-06-01
- Migrate to Vidarr

## 1.0.3 - 2021-11-04
- Workflow input is now fastq pairs, if multiple fastq pairs are defined the workflow will align each pair to the reference genome using bwaMem and then perform a samtools merge. The next steps remain the same. 
