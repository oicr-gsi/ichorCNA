# Changelog
All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased] - 2025-04-23
### Added
[GRD-795](https://jira.oicr.on.ca/browse/GRD-795) - Expanded built-in documentation (metadata changes only).

## [1.3.0] - 2024-06-25
### Added
[GRD-797](https://jira.oicr.on.ca/browse/GRD-797)] - add vidarr labels to outputs (changes to metadata only)

## [1.2.0] - 2023-08-10
### Added
- Assembly-specific modules are specified inside the workflow

## [1.1.2] - 2023-04-12
### Changed
- Update the pdfs provisioned out

## [1.1.1] - 2022-09-15
### Added
- Add bamQC as subworkflow to produce metrics on bam file used for ichorCNA
### Changed
- Fix typos in json metrics file

## [1.1.0] - 2022-07-27
### Changed
- Collect additional metrics
- Output metrics in json file
- Output all plots with annotations

## [1.0.7] - 2022-02-18
### Added
- Add output prefix to bwamem task
### Changed
- Output coverage report always
- Split index and calculate coverage tasks

## [1.0.6] - 2022-01-21
### Added
- Add trimming to alignment step
### Changed
- Update regression test module

## [1.0.5] - 2022-01-07
### Changed
- Workflow input can be fastq pairs or an array of bam files. There is also an option to provision out the bam file used for the analysis, which includes a json file with coverage information.

## [Unreleased] - 2021-11-24
### Changed
[GP-2881](https://jira.oicr.on.ca/browse/GP-2881) make regression test more robust

## [1.0.4] - 2021-11-17
### Changed
- Increment version to bypass Jenkins build error

## [1.0.3] - 2021-11-04
### Changed
- Workflow input is now fastq pairs, if multiple fastq pairs are defined the workflow will align each pair to the reference genome using bwaMem and then perform a samtools merge. The next steps remain the same.

## [1.0.2] - 2021-06-01
### Changed
- Migrate to Vidarr
