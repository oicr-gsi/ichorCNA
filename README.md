# ichorCNA

Workflow for estimating the fraction of tumor in cell-free DNA from sWGS

## Overview

## Dependencies

* [samtools 1.9](http://www.htslib.org/)
* [samtools 1.14](http://www.htslib.org/)
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
`outputFileNamePrefix`|String|Output prefix to prefix output file names with.
`windowSize`|Int|The size of non-overlapping windows.
`minimumMappingQuality`|Int|Mapping quality value below which reads are ignored.
`chromosomesToAnalyze`|String|Chromosomes in the bam reference file.
`provisionBam`|Boolean|Boolean, to provision out bam file and coverage metrics
`inputType`|String|one of either fastq or bam
`bwaMem.runBwaMem_bwaRef`|String|The reference genome to align the sample with by BWA
`bwaMem.runBwaMem_modules`|String|Required environment modules


#### Optional workflow parameters:
Parameter|Value|Default|Description
---|---|---|---
`inputGroups`|Array[InputGroup]?|None|Array of fastq files and their read groups (optional).
`inputBam`|Array[File]?|None|Array of one or multiple bam files (optional).


#### Optional task parameters:
Parameter|Value|Default|Description
---|---|---|---
`bwaMem.adapterTrimmingLog_timeout`|Int|48|Hours before task timeout
`bwaMem.adapterTrimmingLog_jobMemory`|Int|12|Memory allocated indexing job
`bwaMem.indexBam_timeout`|Int|48|Hours before task timeout
`bwaMem.indexBam_modules`|String|"samtools/1.9"|Modules for running indexing job
`bwaMem.indexBam_jobMemory`|Int|12|Memory allocated indexing job
`bwaMem.bamMerge_timeout`|Int|72|Hours before task timeout
`bwaMem.bamMerge_modules`|String|"samtools/1.9"|Required environment modules
`bwaMem.bamMerge_jobMemory`|Int|32|Memory allocated indexing job
`bwaMem.runBwaMem_timeout`|Int|96|Hours before task timeout
`bwaMem.runBwaMem_jobMemory`|Int|32|Memory allocated for this job
`bwaMem.runBwaMem_threads`|Int|8|Requested CPU threads
`bwaMem.runBwaMem_addParam`|String?|None|Additional BWA parameters
`bwaMem.adapterTrimming_timeout`|Int|48|Hours before task timeout
`bwaMem.adapterTrimming_jobMemory`|Int|16|Memory allocated for this job
`bwaMem.adapterTrimming_addParam`|String?|None|Additional cutadapt parameters
`bwaMem.adapterTrimming_modules`|String|"cutadapt/1.8.3"|Required environment modules
`bwaMem.slicerR2_timeout`|Int|48|Hours before task timeout
`bwaMem.slicerR2_jobMemory`|Int|16|Memory allocated for this job
`bwaMem.slicerR2_modules`|String|"slicer/0.3.0"|Required environment modules
`bwaMem.slicerR1_timeout`|Int|48|Hours before task timeout
`bwaMem.slicerR1_jobMemory`|Int|16|Memory allocated for this job
`bwaMem.slicerR1_modules`|String|"slicer/0.3.0"|Required environment modules
`bwaMem.countChunkSize_timeout`|Int|48|Hours before task timeout
`bwaMem.countChunkSize_jobMemory`|Int|16|Memory allocated for this job
`bwaMem.outputFileNamePrefix`|String|"output"|Prefix for output file
`bwaMem.numChunk`|Int|1|number of chunks to split fastq file [1, no splitting]
`bwaMem.doTrim`|Boolean|false|if true, adapters will be trimmed before alignment
`bwaMem.trimMinLength`|Int|1|minimum length of reads to keep [1]
`bwaMem.trimMinQuality`|Int|0|minimum quality of read ends to keep [0]
`bwaMem.adapter1`|String|"AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC"|adapter sequence to trim from read 1 [AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC]
`bwaMem.adapter2`|String|"AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT"|adapter sequence to trim from read 2 [AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT]
`bamMerge.jobMemory`|Int|32|Memory allocated indexing job
`bamMerge.modules`|String|"samtools/1.9"|Required environment modules
`bamMerge.timeout`|Int|72|Hours before task timeout
`inputBamMerge.jobMemory`|Int|32|Memory allocated indexing job
`inputBamMerge.modules`|String|"samtools/1.9"|Required environment modules
`inputBamMerge.timeout`|Int|72|Hours before task timeout
`calculateCoverage.jobMemory`|Int|8|Memory (in GB) to allocate to the job.
`calculateCoverage.modules`|String|"samtools/1.14"|Environment module name and version to load (space separated) before command execution.
`calculateCoverage.timeout`|Int|12|Maximum amount of time (in hours) the task can run for.
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
`bam`|File?|alignment file in bam format used for the analysis (merged if input is multiple fastqs or bams).
`bamIndex`|File?|output index file for bam aligned to genome.
`coverageReport`|File?|json file with the mean coverage for outbam.
`segments`|File|Segments called by the Viterbi algorithm.  Format is compatible with IGV.
`segmentsWithSubclonalStatus`|File|Same as `segments` but also includes subclonal status of segments (0=clonal, 1=subclonal). Format not compatible with IGV.
`estimatedCopyNumber`|File|Estimated copy number, log ratio, and subclone status for each bin/window.
`convergedParameters`|File|Final converged parameters for optimal solution. Also contains table of converged parameters for all solutions.
`correctedDepth`|File|Log2 ratio of each bin/window after correction for GC and mappability biases.
`rData`|File|Saved R image after ichorCNA has finished. Results for all solutions will be included.
`plots`|File|Archived directory of plots.


## Commands
 This section lists command(s) run by ichorCNA workflow

 * Running ichorCNA workflow

 IchorCNA allows for quantification of tumor content in cfDNA. The input for this workflow is an array of fastq pairs with their read group information or an array of bam file(s). When the input is fastqs, the ichorCNA workflow first calls bwaMem for an alignment to the specified reference genome; then if multiple fastq pairs are specified the bam files are merged using samtools. When the input is bam(s) the workflow will merge if necessary or continue. The next step prepares the data for ichorCNA which is the final step in the workflow.

 MERGE BAMS
  ```
 samtools merge -c ~{resultMergedBam} ~{sep=" " bams}
  ```
 CALCULATE COVERAGE
 ```
 samtools index ~{inputbam} ~{resultBai}

 samtools coverage ~{inputbam} | grep -P "^chr\d+\t|^chrX\t|^chrY\t" | awk '{ space += ($3-$2)+1; bases += $7*($3-$2);} END { print bases/space }' | awk '{print "{\"mean coverage\":" $1 "}"}' > ~{outputFileNamePrefix}_coverage.json
 ```
 READCOUNTER
  ```
 samtools index ~{bam}

 # calculate chromosomes to analyze (with reads) from input data
 CHROMOSOMES_WITH_READS=$(samtools view ~{bam} $(tr ',' ' ' <<< ~{chromosomesToAnalyze}) | cut -f3 | sort -V | uniq | paste -s -d, -)

 # write out a chromosomes with reads for ichorCNA
 # split onto new lines (for wdl read_lines), exclude chrY, remove chr prefix, wrap in single quotes for ichorCNA
 echo "${CHROMOSOMES_WITH_READS}" | tr ',' '\n' | grep -v chrY | sed "s/chr//g" | sed -e "s/\(.*\)/'\1'/" > ichorCNAchrs.txt

 # convert
 readCounter \
 --window ~{windowSize} \
 --quality ~{minimumMappingQuality} \
 --chromosome "${CHROMOSOMES_WITH_READS}" \
 ~{bam} | sed "s/chrom=chr/chrom=/" > ~{outputFileNamePrefix}.wig
  ```
 RUN ICHORCNA
  ```
     runIchorCNA \
     --WIG ~{wig} \
     ~{"--NORMWIG " + normalWig} \
     --gcWig ~{gcWig} \
     ~{"--mapWig " + mapWig} \
     ~{"--normalPanel " + normalPanel} \
     ~{"--exons.bed " + exonsBed} \
     --id ~{outputFileNamePrefix} \
     ~{"--centromere " + centromere} \
     ~{"--minMapScore " + minMapScore} \
     ~{"--rmCentromereFlankLength " + rmCentromereFlankLength} \
     ~{"--normal " + normal} \
     ~{"--scStates " + scStates} \
     ~{"--coverage " + coverage} \
     ~{"--lambda " + lambda} \
     ~{"--lambdaScaleHyperParam " + lambdaScaleHyperParam} \
     ~{"--ploidy " + ploidy} \
     ~{"--maxCN " + maxCN} \
     ~{true="--estimateNormal True" false="--estimateNormal False" estimateNormal} \
     ~{true="--estimateScPrevalence True" false="--estimateScPrevalence  False" estimateScPrevalence} \
     ~{true="--estimatePloidy True" false="--estimatePloidy False" estimatePloidy} \
     ~{"--maxFracCNASubclone " + maxFracCNASubclone} \
     ~{"--maxFracGenomeSubclone " + maxFracGenomeSubclone} \
     ~{"--minSegmentBins " + minSegmentBins} \
     ~{"--altFracThreshold " + altFracThreshold} \
     ~{"--chrNormalize " + chrNormalize} \
     ~{"--chrTrain " + chrTrain} \
     --chrs "c(~{sep="," chrs})" \
     ~{"--genomeBuild " + genomeBuild} \
     ~{"--genomeStyle " + genomeStyle} \
     ~{true="--normalizeMaleX True" false="--normalizeMaleX False" normalizeMaleX} \
     ~{"--fracReadsInChrYForMale " + fracReadsInChrYForMale} \
     ~{true="--includeHOMD True" false="--includeHOMD False" includeHOMD} \
     ~{"--txnE " + txnE} \
     ~{"--txnStrength " + txnStrength} \
     ~{"--plotFileType " + plotFileType} \
     ~{"--plotYLim " + plotYLim} \
     ~{"--libdir " + libdir} \
     --outDir ~{outDir}

     # compress directory of plots
     tar -zcvf "~{outputFileNamePrefix}_plots.tar.gz" "~{outputFileNamePrefix}"
  ```
 ## Support

For support, please file an issue on the [Github project](https://github.com/oicr-gsi) or send an email to gsi@oicr.on.ca .

_Generated with generate-markdown-readme (https://github.com/oicr-gsi/gsi-wdl-tools/)_
