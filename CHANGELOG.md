## 1.2.0 - 2023-08-10
- Assembly-specific modules are specified inside the workflow

## 1.1.2 - 2023-04-12
- Update the pdfs provisioned out

## 1.1.1 - 2022-09-15
- Fix typos in json metrics file
- Add bamQC as subworkflow to produce metrics on bam file used for ichorCNA

## 1.1.0 - 2022-07-27
- Collect additional metrics
- Output metrics in json file
- Output all plots with annotations

## 1.0.7 - 2022-02-18
- Output coverage report always
- Split index and calculate coverage tasks
- Add output prefix to bwamem task

## 1.0.6 - 2022-01-21
- Add trimming to alignment step
- Update regression test module

## 1.0.5 - 2022-01-07
- Workflow input can be fastq pairs or an array of bam files. There is also an option to provision out the bam file used for the analysis, which includes a json file with coverage information.

## Unreleased - 2021-11-24
[GP-2881](https://jira.oicr.on.ca/browse/GP-2881) make regression test more robust

## 1.0.4 - 2021-11-17
- Increment version to bypass Jenkins build error

## 1.0.3 - 2021-11-04
- Workflow input is now fastq pairs, if multiple fastq pairs are defined the workflow will align each pair to the reference genome using bwaMem and then perform a samtools merge. The next steps remain the same.

## 1.0.2 - 2021-06-01
- Migrate to Vidarr
