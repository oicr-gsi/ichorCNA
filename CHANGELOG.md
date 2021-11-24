## 1.0.2 - 2021-06-01
- Migrate to Vidarr

## 1.0.3 - 2021-11-04
- Workflow input is now fastq pairs, if multiple fastq pairs are defined the workflow will align each pair to the reference genome using bwaMem and then perform a samtools merge. The next steps remain the same.

## 1.0.4 - 2021-11-17
- Increment version to bypass Jenkins build error

## Unreleased - 2021-11-24
[GP-2881](https://jira.oicr.on.ca/browse/GP-2881) make regression test more robust 
