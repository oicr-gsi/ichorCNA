# ichorCNA

Workflow for estimating the fraction of tumor in cell-free DNA from sWGS

## Overview

## Dependencies

* [samtools 1.9](http://www.htslib.org/)
* [samtools 1.14](http://www.htslib.org/)
* [hmmcopy-utils 0.1.1](https://github.com/broadinstitute/ichorCNA)
* [ichorcna 0.2](https://shahlab.ca/projects/hmmcopy_utils/)
* [python 3.9](python.org)
* [pandas 1.4.2](https://pandas.pydata.org/)


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
`reference`|String|Genome build (e.g. "hg19" or "hg38")
`bwaMem.runBwaMem_bwaRef`|String|The reference genome to align the sample with by BWA
`bwaMem.runBwaMem_modules`|String|Required environment modules
`bamQC.bamQCMetrics_workflowVersion`|String|Workflow version string
`bamQC.bamQCMetrics_refSizesBed`|String|Path to human genome BED reference with chromosome sizes
`bamQC.bamQCMetrics_refFasta`|String|Path to human genome FASTA reference
`bamQC.metadata`|Map[String,String]|JSON file containing metadata


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
`bwaMem.numChunk`|Int|1|number of chunks to split fastq file [1, no splitting]
`bwaMem.trimMinLength`|Int|1|minimum length of reads to keep [1]
`bwaMem.trimMinQuality`|Int|0|minimum quality of read ends to keep [0]
`bwaMem.adapter1`|String|"AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC"|adapter sequence to trim from read 1 [AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC]
`bwaMem.adapter2`|String|"AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT"|adapter sequence to trim from read 2 [AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT]
`preMergeBamMetricsFastqInput.jobMemory`|Int|8|Memory (in GB) to allocate to the job.
`preMergeBamMetricsFastqInput.modules`|String|"samtools/1.14"|Environment module name and version to load (space separated) before command execution.
`preMergeBamMetricsFastqInput.timeout`|Int|12|Maximum amount of time (in hours) the task can run for.
`bamMerge.jobMemory`|Int|32|Memory allocated indexing job
`bamMerge.modules`|String|"samtools/1.9"|Required environment modules
`bamMerge.timeout`|Int|72|Hours before task timeout
`preMergeBamMetrics.jobMemory`|Int|8|Memory (in GB) to allocate to the job.
`preMergeBamMetrics.modules`|String|"samtools/1.14"|Environment module name and version to load (space separated) before command execution.
`preMergeBamMetrics.timeout`|Int|12|Maximum amount of time (in hours) the task can run for.
`inputBamMerge.jobMemory`|Int|32|Memory allocated indexing job
`inputBamMerge.modules`|String|"samtools/1.9"|Required environment modules
`inputBamMerge.timeout`|Int|72|Hours before task timeout
`indexBam.jobMemory`|Int|8|Memory (in GB) to allocate to the job.
`indexBam.modules`|String|"samtools/1.9"|Environment module name and version to load (space separated) before command execution.
`indexBam.timeout`|Int|12|Maximum amount of time (in hours) the task can run for.
`runReadCounter.mem`|Int|8|Memory (in GB) to allocate to the job.
`runReadCounter.modules`|String|"samtools/1.9 hmmcopy-utils/0.1.1"|Environment module name and version to load (space separated) before command execution.
`runReadCounter.timeout`|Int|12|Maximum amount of time (in hours) the task can run for.
`runIchorCNA.normalWig`|File?|None|Normal WIG file. Default: [NULL].
`runIchorCNA.gcWig`|String|None|GC-content WIG file.
`runIchorCNA.mapWig`|String|None|Mappability score WIG file. Default: [NULL].
`runIchorCNA.normalPanel`|String|None|Median corrected depth from panel of normals. Default: [NULL].
`runIchorCNA.exonsBed`|String?|None|Bed file containing exon regions. Default: [NULL].
`runIchorCNA.centromere`|String|None|File containing Centromere locations; if not provided then will use hg19 version from ichorCNA package.
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
`runIchorCNA.genomeBuild`|String|None|Genome build.
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
`bamQC.collateResults_timeout`|Int|1|hours before task timeout
`bamQC.collateResults_threads`|Int|4|Requested CPU threads
`bamQC.collateResults_jobMemory`|Int|8|Memory allocated for this job
`bamQC.collateResults_modules`|String|"python/3.6"|required environment modules
`bamQC.cumulativeDistToHistogram_timeout`|Int|1|hours before task timeout
`bamQC.cumulativeDistToHistogram_threads`|Int|4|Requested CPU threads
`bamQC.cumulativeDistToHistogram_jobMemory`|Int|8|Memory allocated for this job
`bamQC.cumulativeDistToHistogram_modules`|String|"python/3.6"|required environment modules
`bamQC.runMosdepth_timeout`|Int|4|hours before task timeout
`bamQC.runMosdepth_threads`|Int|4|Requested CPU threads
`bamQC.runMosdepth_jobMemory`|Int|16|Memory allocated for this job
`bamQC.runMosdepth_modules`|String|"mosdepth/0.2.9"|required environment modules
`bamQC.bamQCMetrics_timeout`|Int|4|hours before task timeout
`bamQC.bamQCMetrics_threads`|Int|4|Requested CPU threads
`bamQC.bamQCMetrics_jobMemory`|Int|16|Memory allocated for this job
`bamQC.bamQCMetrics_modules`|String|"bam-qc-metrics/0.2.5"|required environment modules
`bamQC.bamQCMetrics_normalInsertMax`|Int|1500|Maximum of expected insert size range
`bamQC.markDuplicates_timeout`|Int|4|hours before task timeout
`bamQC.markDuplicates_threads`|Int|4|Requested CPU threads
`bamQC.markDuplicates_jobMemory`|Int|16|Memory allocated for this job
`bamQC.markDuplicates_modules`|String|"picard/2.21.2"|required environment modules
`bamQC.markDuplicates_picardMaxMemMb`|Int|6000|Memory requirement in MB for running Picard JAR
`bamQC.markDuplicates_opticalDuplicatePixelDistance`|Int|100|Maximum offset between optical duplicate clusters
`bamQC.downsampleRegion_timeout`|Int|4|hours before task timeout
`bamQC.downsampleRegion_threads`|Int|4|Requested CPU threads
`bamQC.downsampleRegion_jobMemory`|Int|16|Memory allocated for this job
`bamQC.downsampleRegion_modules`|String|"samtools/1.9"|required environment modules
`bamQC.downsample_timeout`|Int|4|hours before task timeout
`bamQC.downsample_threads`|Int|4|Requested CPU threads
`bamQC.downsample_jobMemory`|Int|16|Memory allocated for this job
`bamQC.downsample_modules`|String|"samtools/1.9"|required environment modules
`bamQC.downsample_randomSeed`|Int|42|Random seed for pre-downsampling (if any)
`bamQC.downsample_downsampleSuffix`|String|"downsampled.bam"|Suffix for output file
`bamQC.findDownsampleParamsMarkDup_timeout`|Int|4|hours before task timeout
`bamQC.findDownsampleParamsMarkDup_threads`|Int|4|Requested CPU threads
`bamQC.findDownsampleParamsMarkDup_jobMemory`|Int|16|Memory allocated for this job
`bamQC.findDownsampleParamsMarkDup_modules`|String|"python/3.6"|required environment modules
`bamQC.findDownsampleParamsMarkDup_customRegions`|String|""|Custom downsample regions; overrides chromosome and interval parameters
`bamQC.findDownsampleParamsMarkDup_intervalStart`|Int|100000|Start of interval in each chromosome, for very large BAMs
`bamQC.findDownsampleParamsMarkDup_baseInterval`|Int|15000|Base width of interval in each chromosome, for very large BAMs
`bamQC.findDownsampleParamsMarkDup_chromosomes`|Array[String]|["chr12", "chr13", "chrXII", "chrXIII"]|Array of chromosome identifiers for downsampled subset
`bamQC.findDownsampleParamsMarkDup_threshold`|Int|10000000|Minimum number of reads to conduct downsampling
`bamQC.findDownsampleParams_timeout`|Int|4|hours before task timeout
`bamQC.findDownsampleParams_threads`|Int|4|Requested CPU threads
`bamQC.findDownsampleParams_jobMemory`|Int|16|Memory allocated for this job
`bamQC.findDownsampleParams_modules`|String|"python/3.6"|required environment modules
`bamQC.findDownsampleParams_preDSMultiplier`|Float|1.5|Determines target size for pre-downsampled set (if any). Must have (preDSMultiplier) < (minReadsRelative).
`bamQC.findDownsampleParams_precision`|Int|8|Number of decimal places in fraction for pre-downsampling
`bamQC.findDownsampleParams_minReadsRelative`|Int|2|Minimum value of (inputReads)/(targetReads) to allow pre-downsampling
`bamQC.findDownsampleParams_minReadsAbsolute`|Int|10000|Minimum value of targetReads to allow pre-downsampling
`bamQC.findDownsampleParams_targetReads`|Int|100000|Desired number of reads in downsampled output
`bamQC.indexBamFile_timeout`|Int|4|hours before task timeout
`bamQC.indexBamFile_threads`|Int|4|Requested CPU threads
`bamQC.indexBamFile_jobMemory`|Int|16|Memory allocated for this job
`bamQC.indexBamFile_modules`|String|"samtools/1.9"|required environment modules
`bamQC.countInputReads_timeout`|Int|4|hours before task timeout
`bamQC.countInputReads_threads`|Int|4|Requested CPU threads
`bamQC.countInputReads_jobMemory`|Int|16|Memory allocated for this job
`bamQC.countInputReads_modules`|String|"samtools/1.9"|required environment modules
`bamQC.updateMetadata_timeout`|Int|4|hours before task timeout
`bamQC.updateMetadata_threads`|Int|4|Requested CPU threads
`bamQC.updateMetadata_jobMemory`|Int|16|Memory allocated for this job
`bamQC.updateMetadata_modules`|String|"python/3.6"|required environment modules
`bamQC.filter_timeout`|Int|4|hours before task timeout
`bamQC.filter_threads`|Int|4|Requested CPU threads
`bamQC.filter_jobMemory`|Int|16|Memory allocated for this job
`bamQC.filter_modules`|String|"samtools/1.9"|required environment modules
`bamQC.filter_minQuality`|Int|30|Minimum alignment quality to pass filter
`getMetrics.jobMemory`|Int|8|Memory (in GB) to allocate to the job.
`getMetrics.modules`|String|"samtools/1.14"|Environment module name and version to load (space separated) before command execution.
`getMetrics.timeout`|Int|12|Maximum amount of time (in hours) the task can run for.
`createJson.jobMemory`|Int|8|Memory (in GB) to allocate to the job.
`createJson.modules`|String|"pandas/1.4.2"|Environment module name and version to load (space separated) before command execution.
`createJson.timeout`|Int|12|Maximum amount of time (in hours) the task can run for.


### Outputs

Output | Type | Description
---|---|---
`genomeWideAll`|Pair[File,Map[String,String]]|Genome wide plots for each solution
`genomeWide`|Pair[File,Map[String,String]]|Genome wide plots for the selected solution
`bam`|File?|Bam file used as input to ichorCNA (only produced when provisionBam is True)
`bamIndex`|File?|Bam index for bam file used as input to ichorCNA (only produced when provisionBam is True)
`jsonMetrics`|File|Report on bam coverage, read counts and ichorCNA metrics.
`segments`|File|Segments called by the Viterbi algorithm.  Format is compatible with IGV.
`segmentsWithSubclonalStatus`|File|Same as segments but also includes subclonal status of segments (0=clonal, 1=subclonal). Format not compatible with IGV.
`estimatedCopyNumber`|File|Estimated copy number, log ratio, and subclone status for each bin/window.
`convergedParameters`|File|Final converged parameters for optimal solution. Also contains table of converged parameters for all solutions.
`correctedDepth`|File|Log2 ratio of each bin/window after correction for GC and mappability biases.
`rData`|File|Saved R image after ichorCNA has finished. Results for all solutions will be included.
`plots`|File|Archived directory of plots.
`bamQCresult`|File|bamQC report.


## Commands
 This section lists command(s) run by ichorCNA workflow
 
 * Running ichorCNA workflow
 
 IchorCNA allows for quantification of tumor content in cfDNA. The input for this workflow is an array of fastq pairs with their read group information. This ichorCNA workflow first calls bwaMem for an alignment to the specified reference genome; then if multiple fastq pairs are specified the bam files are merged using samtools. The next step prepares the data for ichorCNA which is the final step in the workflow.
 
 MERGE BAMS
 ```
 samtools merge \
 -c \
 ~{resultMergedBam} \
 ~{sep=" " bams}
 ```
 COLLECT PRE-MERGE BAM METRICS
 ```
 echo run,read_count > ~{outputFileNamePrefix}_pre_merge_bam_metrics.csv
 for file in ~{sep=' ' bam}
 do
   run=$(samtools view -H "${file}" | grep '^@RG' | cut -f 2 | cut -f 2 -d ":" | cut -f 1 -d "-")
   read_count=$(samtools stats "${file}" | grep ^SN | grep "raw total sequences" | cut -f 3)
   echo $run,$read_count >> ~{outputFileNamePrefix}_pre_merge_bam_metrics.csv
 done;
 ``` 
 INDEX BAM
 ```
 samtools index ~{inputbam} ~{resultBai}
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
 
     #create txt file with plot full path
     ls $PWD/~{outputFileNamePrefix}/*genomeWide_n* > "~{outputFileNamePrefix}"_plots.txt
  ```
  COLLECT FINAL BAM AND ICHORCNA METRICS
  ```
  echo coverage,read_count,tumor_fraction,ploidy > ~{outputFileNamePrefix}_bam_metrics.csv
  coverage=$(samtools coverage ~{inputbam} | grep -P "^chr\d+\t|^chrX\t|^chrY\t" | awk '{ space += ($3-$2)+1; bases += $7*($3-$2);} END { print bases/space }')
  read_count=$(samtools stats ~{inputbam} | head -n 8 | tail -n 1 | cut -f 3)
  tumor_fraction=$(cat ~{params} | head -n 2 | tail -n 1 | cut -f 2)
  ploidy=$(cat ~{params} | head -n 2 | tail -n 1 | cut -f 3)
  echo $coverage,$read_count,$tumor_fraction,$ploidy >> ~{outputFileNamePrefix}_bam_metrics.csv
  cat ~{params} | tail -n 17 > ~{outputFileNamePrefix}_all_sols_metrics.csv
  ```
 CREATE JSON WITH METRICS COLLECTED
 ```
 python3 <<CODE
 import csv, json
 import pandas as pd
 
 ### create json file with all metrics
 
 bam_metric = pd.read_csv("~{bamMetrics}")
 pre_metric = pd.read_csv("~{preBamMetrics}")
 all_sols = pd.read_csv("~{allSolsMetrics}", sep="\t")
 all_sols["tumor_fraction"] = round(1 - all_sols["n_est"],3)
 all_sols["solution"] = all_sols["init"]
 pre_metric_dict = pre_metric.to_dict('index')
 bam_metric_dict = bam_metric.to_dict('records')[0]
 with open("~{plotsFile}") as f:
   lines = f.readlines()
 
 #reorganize lane sequencing data
 lanes = []
 for lane in pre_metric_dict:
   lanes.append(pre_metric_dict[lane])
 
 #find selected solution
 selected_sol = ""
 for index, row in all_sols.iterrows():
   if round(row["tumor_fraction"],2) == round(bam_metric_dict["tumor_fraction"],2) and row["phi_est"] == bam_metric_dict["ploidy"]:
     selected_sol = row["init"]
 
 #selecting metrics from all solutions
 all_sols_metrics = {}
 for index, row in all_sols.iterrows():
   all_sols_metrics[row["solution"]] = {"tumor_fraction":row["tumor_fraction"],
                                       "ploidy":row["phi_est"],
                                       "loglik":row["loglik"]}
 
 metrics_dict = {"mean_coverage": bam_metric_dict["coverage"],
                 "total_reads": bam_metric_dict["read_count"],
                 "lanes_sequenced": len(pre_metric_dict),
                 "reads_per_lane":lanes,
                 "best_solution": selected_sol,
                 "tumor_fraction": bam_metric_dict["tumor_fraction"],
                 "ploidy": bam_metric_dict["ploidy"],
                 "solutions": all_sols_metrics}
 
 with open("~{outputFileNamePrefix}_metrics.json", "w") as outfile:
   json.dump(metrics_dict, outfile)
 
 ### create json output file for annotations
 output_list = []
 for line in lines:
    pdf_dict = {}
    line = line.strip()
    pdf_dict["left"] = line
    pdf_dict["right"] = {}
    pdf_dict["right"]["tumor_fraction"] = bam_metric_dict["tumor_fraction"]
    pdf_dict["right"]["ploidy"] = bam_metric_dict["ploidy"]
    output_list.append(pdf_dict)
 output_dict = {}
 output_dict["pdfs"] = output_list
 
 with open("~{outputFileNamePrefix}_outputs.json", "w") as outPdfJson:
     json.dump(output_dict, outPdfJson)
 
 CODE
 ```
 ## Support

For support, please file an issue on the [Github project](https://github.com/oicr-gsi) or send an email to gsi@oicr.on.ca .

_Generated with generate-markdown-readme (https://github.com/oicr-gsi/gsi-wdl-tools/)_
