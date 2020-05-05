# ichorCNA

Workflow for estimating the fraction of tumor in cell-free DNA from sWGS

## Overview

## Dependencies

* [samtools 1.9](http://www.htslib.org/)
* [hmmcopy-utils 0.1.1](https://github.com/broadinstitute/ichorCNA)
* [ichorcna 0.2](https://shahlab.ca/projects/hmmcopy_utils/)


## Usage

### Cromwell
```
java -jar cromwell.jar run ichorCNA.wdl --inputs inputs.json
```

### Inputs

#### Required workflow parameters:
Parameter|Value|Description
---|---|---
`bam`|File|Input bam.
`bamIndex`|File|Input bam index (must be .bam.bai).
`windowSize`|Int|The size of non-overlapping windows.
`minimumMappingQuality`|Int|Mapping quality value below which reads are ignored.
`chromosomesToAnalyze`|String|Chromosomes in the bam reference file.


#### Optional workflow parameters:
Parameter|Value|Default|Description
---|---|---|---
`outputFileNamePrefix`|String|basename(bam,'.bam')|Output prefix to prefix output file names with.


#### Optional task parameters:
Parameter|Value|Default|Description
---|---|---|---
`runReadCounter.mem`|Int|8|Memory (in GB) to allocate to the job.
`runReadCounter.modules`|String|"samtools/1.9 hmmcopy-utils/0.1.1"|Environment module name and version to load (space separated) before command execution.
`runReadCounter.timeout`|Int|12|Maximum amount of time (in hours) the task can run for.
`runIchorCNA.normalWig`|File?|None|Normal WIG file. Default: [NULL].
`runIchorCNA.gcWig`|String|"$ICHORCNA_ROOT/lib/R/ichorCNA/extdata/gc_hg19_1000kb.wig"|GC-content WIG file.
`runIchorCNA.mapWig`|String|"$ICHORCNA_ROOT/lib/R/ichorCNA/extdata/map_hg19_1000kb.wig"|Mappability score WIG file. Default: [NULL].
`runIchorCNA.normalPanel`|String|"$ICHORCNA_ROOT/lib/R/ichorCNA/extdata/HD_ULP_PoN_1Mb_median_normAutosome_mapScoreFiltered_median.rds"|Median corrected depth from panel of normals. Default: [NULL].
`runIchorCNA.exonsBed`|String?|None|Bed file containing exon regions. Default: [NULL].
`runIchorCNA.centromere`|String|"$ICHORCNA_ROOT/lib/R/ichorCNA/extdata/GRCh37.p13_centromere_UCSC-gapTable.txt"|File containing Centromere locations; if not provided then will use hg19 version from ichorCNA package.
`runIchorCNA.minMapScore`|Float?|None|Include bins with a minimum mappability score of this value. Default: [0.9].
`runIchorCNA.rmCentromereFlankLength`|Int?|None|Length of region flanking centromere to remove. Default: [1e+05].
`runIchorCNA.normal`|String|"\"c(0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9)\""|Initial normal contamination; can be more than one value if additional normal initializations are desired. Default: [0.5]
`runIchorCNA.scStates`|String|"\"c(1, 3)\""|Subclonal states to consider.
`runIchorCNA.coverage`|String?|None|PICARD sequencing coverage.
`runIchorCNA.lambda`|String?|None|Initial Student's t precision; must contain 4 values (e.g. c(1500,1500,1500,1500)); if not provided then will automatically use based on variance of data.
`runIchorCNA.lambdaScaleHyperParam`|Int?|None|Hyperparameter (scale) for Gamma prior on Student's-t precision. Default: [3].
`runIchorCNA.ploidy`|String|"\"c(2,3)\""|Initial tumour ploidy; can be more than one value if additional ploidy initializations are desired. Default: [2]
`runIchorCNA.maxCN`|Int|5|Total clonal CN states.
`runIchorCNA.estimateNormal`|Boolean|true|Estimate normal?
`runIchorCNA.estimateScPrevalence`|Boolean|true|Estimate subclonal prevalence?
`runIchorCNA.estimatePloidy`|Boolean|true|Estimate tumour ploidy?
`runIchorCNA.maxFracCNASubclone`|Float?|None|Exclude solutions with fraction of subclonal events greater than this value. Default: [0.7].
`runIchorCNA.maxFracGenomeSubclone`|Float?|None|Exclude solutions with subclonal genome fraction greater than this value. Default: [0.5].
`runIchorCNA.minSegmentBins`|String?|None|Minimum number of bins for largest segment threshold required to estimate tumor fraction; if below this threshold, then will be assigned zero tumor fraction.
`runIchorCNA.altFracThreshold`|Float?|None|Minimum proportion of bins altered required to estimate tumor fraction; if below this threshold, then will be assigned zero tumor fraction. Default: [0.05].
`runIchorCNA.chrNormalize`|String?|None|Specify chromosomes to normalize GC/mappability biases. Default: [c(1:22)].
`runIchorCNA.chrTrain`|String|"\"c(1:22)\""|Specify chromosomes to estimate params. Default: [c(1:22)].
`runIchorCNA.genomeBuild`|String?|None|Genome build. Default: [hg19].
`runIchorCNA.genomeStyle`|String?|None|NCBI or UCSC chromosome naming convention; use UCSC if desired output is to have "chr" string. [Default: NCBI].
`runIchorCNA.normalizeMaleX`|Boolean?|None|If male, then normalize chrX by median. Default: [TRUE].
`runIchorCNA.fracReadsInChrYForMale`|Float?|None|Threshold for fraction of reads in chrY to assign as male. Default: [0.001].
`runIchorCNA.includeHOMD`|Boolean|true|If FALSE, then exclude HOMD state. Useful when using large bins (e.g. 1Mb). Default: [FALSE].
`runIchorCNA.txnE`|Float|0.9999|Self-transition probability. Increase to decrease number of segments. Default: [0.9999999]
`runIchorCNA.txnStrength`|Int|10000|Transition pseudo-counts. Exponent should be the same as the number of decimal places of --txnE. Default: [1e+07].
`runIchorCNA.plotFileType`|String?|None|File format for output plots. Default: [pdf].
`runIchorCNA.plotYLim`|String?|None|ylim to use for chromosome plots. Default: [c(-2,2)].
`runIchorCNA.outDir`|String|"./"|Output Directory. Default: [./].
`runIchorCNA.libdir`|String?|None|Script library path.
`runIchorCNA.modules`|String|"ichorcna/0.2"|Environment module name and version to load (space separated) before command execution.
`runIchorCNA.mem`|Int|8|Memory (in GB) to allocate to the job.
`runIchorCNA.timeout`|Int|12|Maximum amount of time (in hours) the task can run for.


### Outputs

Output | Type | Description
---|---|---
`segments`|File|Segments called by the Viterbi algorithm.  Format is compatible with IGV.
`segmentsWithSubclonalStatus`|File|Same as `segments` but also includes subclonal status of segments (0=clonal, 1=subclonal). Format not compatible with IGV.
`estimatedCopyNumber`|File|Estimated copy number, log ratio, and subclone status for each bin/window.
`convergedParameters`|File|Final converged parameters for optimal solution. Also contains table of converged parameters for all solutions.
`correctedDepth`|File|Log2 ratio of each bin/window after correction for GC and mappability biases.
`rData`|File|Saved R image after ichorCNA has finished. Results for all solutions will be included.
`plots`|File|Archived directory of plots.


## Niassa + Cromwell

This WDL workflow is wrapped in a Niassa workflow (https://github.com/oicr-gsi/pipedev/tree/master/pipedev-niassa-cromwell-workflow) so that it can used with the Niassa metadata tracking system (https://github.com/oicr-gsi/niassa).

* Building
```
mvn clean install
```

* Testing
```
mvn clean verify \
-Djava_opts="-Xmx1g -XX:+UseG1GC -XX:+UseStringDeduplication" \
-DrunTestThreads=2 \
-DskipITs=false \
-DskipRunITs=false \
-DworkingDirectory=/path/to/tmp/ \
-DschedulingHost=niassa_oozie_host \
-DwebserviceUrl=http://niassa-url:8080 \
-DwebserviceUser=niassa_user \
-DwebservicePassword=niassa_user_password \
-Dcromwell-host=http://cromwell-url:8000
```

## Support

For support, please file an issue on the [Github project](https://github.com/oicr-gsi) or send an email to gsi@oicr.on.ca .

_Generated with generate-markdown-readme (https://github.com/oicr-gsi/gsi-wdl-tools/)_
