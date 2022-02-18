version 1.0

import "imports/pull_bwa.wdl" as bwaMem

struct InputGroup {
  File fastqR1
  File fastqR2
  String readGroups
}

workflow ichorCNA {
  input {
    Array[InputGroup]? inputGroups
    Array[File]? inputBam
    String outputFileNamePrefix
    Int windowSize
    Int minimumMappingQuality
    String chromosomesToAnalyze
    Boolean provisionBam
    String inputType
  }

  parameter_meta {
    inputGroups: "Array of fastq files and their read groups (optional)."
    inputBam: "Array of one or multiple bam files (optional)."
    outputFileNamePrefix: "Output prefix to prefix output file names with."
    windowSize: "The size of non-overlapping windows."
    minimumMappingQuality: "Mapping quality value below which reads are ignored."
    chromosomesToAnalyze: "Chromosomes in the bam reference file."
    provisionBam: "Boolean, to provision out bam file and coverage metrics"
    inputType: "one of either fastq or bam"
  }

  if(inputType=="fastq" && defined(inputGroups)){
    Array[InputGroup] inputGroups_ = select_first([inputGroups])
    scatter (ig in inputGroups_) {
      call bwaMem.bwaMem {
        input:
          fastqR1 = ig.fastqR1,
          fastqR2 = ig.fastqR2,
          readGroups = ig.readGroups,
          doTrim = true,
          outputFileNamePrefix = outputFileNamePrefix
      }
    }

    if (length(inputGroups_) > 1 ) {
      call bamMerge {
        input:
          bams = bwaMem.bwaMemBam,
          outputFileNamePrefix = outputFileNamePrefix
      }
    }

    if (length(inputGroups_) == 1 ) {
      File bwaMemBam = bwaMem.bwaMemBam[0]
    }
  }

  if(inputType=="bam" && defined(inputBam)){
    Array[File] inputBam_ = select_first([inputBam,[]])
    if (length(inputBam_) > 1 ) {
      call bamMerge as inputBamMerge {
        input:
          bams = inputBam_,
          outputFileNamePrefix = outputFileNamePrefix
      }
    }

    if (length(inputBam_) == 1 ) {
      File singleInputBam = inputBam_[0]
    }
  }

  call calculateCoverage {
    input:
      inputbam = select_first([bamMerge.outputMergedBam,bwaMemBam,inputBamMerge.outputMergedBam,singleInputBam]),
      outputFileNamePrefix = outputFileNamePrefix
  }

  if(provisionBam==true){
    call indexBam {
      input:
        inputbam = select_first([bamMerge.outputMergedBam,bwaMemBam,inputBamMerge.outputMergedBam,singleInputBam])
    }
  }

  call runReadCounter{
    input:
      bam= select_first([bamMerge.outputMergedBam,bwaMemBam,inputBamMerge.outputMergedBam,singleInputBam]),
      outputFileNamePrefix=outputFileNamePrefix,
      windowSize=windowSize,
      minimumMappingQuality=minimumMappingQuality,
      chromosomesToAnalyze=chromosomesToAnalyze
  }

  call runIchorCNA {
    input:
      outputFileNamePrefix=outputFileNamePrefix,
      chrs=runReadCounter.ichorCNAchrs,
      wig=runReadCounter.wig
  }

  output {
    File? bam = indexBam.outbam
    File? bamIndex = indexBam.bamIndex
    File coverageReport = calculateCoverage.coverageReport
    File segments = runIchorCNA.segments
    File segmentsWithSubclonalStatus = runIchorCNA.segmentsWithSubclonalStatus
    File estimatedCopyNumber = runIchorCNA.estimatedCopyNumber
    File convergedParameters = runIchorCNA.convergedParameters
    File correctedDepth = runIchorCNA.correctedDepth
    File rData = runIchorCNA.rData
    File plots = runIchorCNA.plots
  }

  meta {
    author: "Michael Laszloffy"
    email: "michael.laszloffy@oicr.on.ca"
    description: "Workflow for estimating the fraction of tumor in cell-free DNA from sWGS"
    dependencies: [
      {
        name: "samtools/1.9",
        url: "http://www.htslib.org/"
      },
      {
        name: "samtools/1.14",
        url: "http://www.htslib.org/"
      },
      {
        name: "hmmcopy-utils/0.1.1",
        url: "https://github.com/broadinstitute/ichorCNA"
      },
      {
        name: "ichorcna/0.2",
        url: "https://shahlab.ca/projects/hmmcopy_utils/"
      }
    ]
  }

}

#copy from bwaMem
task bamMerge{
    input {
        Array[File] bams
        String outputFileNamePrefix
        Int   jobMemory = 32
        String modules  = "samtools/1.9"
        Int timeout     = 72
    }
    parameter_meta {
        bams:  "Input bam files"
        outputFileNamePrefix: "Prefix for output file"
        jobMemory: "Memory allocated indexing job"
        modules:   "Required environment modules"
        timeout:   "Hours before task timeout"
    }

    String resultMergedBam = "~{outputFileNamePrefix}.bam"

    command <<<
        set -euo pipefail
        samtools merge \
        -c \
        ~{resultMergedBam} \
        ~{sep=" " bams}
    >>>

    runtime {
        memory: "~{jobMemory} GB"
        modules: "~{modules}"
        timeout: "~{timeout}"
    }

    output {
        File outputMergedBam = "~{resultMergedBam}"
    }

    meta {
        output_meta: {
            outputMergedBam: "output merged bam aligned to genome"
        }
    }
}

task calculateCoverage {
  input {
    File inputbam
    String outputFileNamePrefix
    Int jobMemory = 8
    String modules = "samtools/1.14"
    Int timeout = 12
  }

  command <<<
  samtools coverage ~{inputbam} | grep -P "^chr\d+\t|^chrX\t|^chrY\t" | awk '{ space += ($3-$2)+1; bases += $7*($3-$2);} END { print bases/space }' | awk '{print "{\"mean coverage\":" $1 "}"}' > ~{outputFileNamePrefix}_coverage.json
  >>>

  output {
  File coverageReport = "~{outputFileNamePrefix}_coverage.json"
  }

  runtime {
    memory: "~{jobMemory} GB"
    modules: "~{modules}"
    timeout: "~{timeout}"
  }

  parameter_meta {
    inputbam: "Input bam."
    outputFileNamePrefix: "Output prefix to prefix output file names with."
    jobMemory: "Memory (in GB) to allocate to the job."
    modules: "Environment module name and version to load (space separated) before command execution."
    timeout: "Maximum amount of time (in hours) the task can run for."
  }

  meta {
    output_meta: {
      coverageReport: "json file with the mean coverage for outbam."
    }
  }
}

task indexBam {
  input {
    File inputbam
    Int jobMemory = 8
    String modules = "samtools/1.9"
    Int timeout = 12
  }

  String resultBai = "~{basename(inputbam)}.bai"

  command <<<
  samtools index ~{inputbam} ~{resultBai}
  >>>

  output {
    File outbam = "~{inputbam}"
    File bamIndex = "~{resultBai}"
  }

  runtime {
    memory: "~{jobMemory} GB"
    modules: "~{modules}"
    timeout: "~{timeout}"
  }

  parameter_meta {
    inputbam: "Input bam."
    jobMemory: "Memory (in GB) to allocate to the job."
    modules: "Environment module name and version to load (space separated) before command execution."
    timeout: "Maximum amount of time (in hours) the task can run for."
  }

  meta {
    output_meta: {
      outbam: "alignment file in bam format used for the analysis (merged if input is multiple fastqs or bams).",
      bamIndex: "output index file for bam aligned to genome.",
    }
  }
}

task runReadCounter {
  input{
    File bam
    String outputFileNamePrefix
    Int windowSize
    Int minimumMappingQuality
    String chromosomesToAnalyze
    Int mem = 8
    String modules = "samtools/1.9 hmmcopy-utils/0.1.1"
    Int timeout = 12
  }

  command <<<
    set -euxo pipefail

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
  >>>

  runtime {
    memory: "~{mem} GB"
    modules: "~{modules}"
    timeout: "~{timeout}"
  }

  output {
    File wig = "~{outputFileNamePrefix}.wig"
    Array[String] ichorCNAchrs = read_lines("ichorCNAchrs.txt")
  }

  parameter_meta {
    bam: "Input bam."
    outputFileNamePrefix: "Output prefix to prefix output file names with."
    windowSize: "The size of non-overlapping windows."
    minimumMappingQuality: "Mapping quality value below which reads are ignored."
    chromosomesToAnalyze: "Chromosomes in the bam reference file."
    mem: "Memory (in GB) to allocate to the job."
    modules: "Environment module name and version to load (space separated) before command execution."
    timeout: "Maximum amount of time (in hours) the task can run for."
  }

  meta {
    output_meta: {
      wig: "Read count file in WIG format",
      ichorCNAchrs: "Chromosomes with reads for ichorCNA (\"chr\" stripped from the name)"
    }
  }
}

task runIchorCNA {
  input {
    String outputFileNamePrefix
    File wig
    File? normalWig
    String gcWig = "$ICHORCNA_ROOT/lib/R/ichorCNA/extdata/gc_hg19_1000kb.wig"
    String mapWig = "$ICHORCNA_ROOT/lib/R/ichorCNA/extdata/map_hg19_1000kb.wig"
    String normalPanel = "$ICHORCNA_ROOT/lib/R/ichorCNA/extdata/HD_ULP_PoN_1Mb_median_normAutosome_mapScoreFiltered_median.rds"
    String? exonsBed
    String centromere = "$ICHORCNA_ROOT/lib/R/ichorCNA/extdata/GRCh37.p13_centromere_UCSC-gapTable.txt"
    Float? minMapScore
    Int? rmCentromereFlankLength
    String normal = "\"c(0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9)\""
    String scStates = "\"c(1, 3)\""
    String? coverage
    String? lambda
    Int? lambdaScaleHyperParam
    String ploidy = "\"c(2,3)\""
    Int maxCN = 5
    Boolean estimateNormal = true
    Boolean estimateScPrevalence = true
    Boolean estimatePloidy = true
    Float? maxFracCNASubclone
    Float? maxFracGenomeSubclone
    String? minSegmentBins
    Float? altFracThreshold
    String? chrNormalize
    String chrTrain = "\"c(1:22)\""
    Array[String] chrs
    String? genomeBuild
    String? genomeStyle
    Boolean? normalizeMaleX
    Float? fracReadsInChrYForMale
    Boolean includeHOMD = true
    Float txnE = 0.9999
    Int txnStrength = 10000
    String? plotFileType
    String? plotYLim
    String outDir = "./"
    String? libdir

    String modules = "ichorcna/0.2"
    Int mem = 8
    Int timeout = 12
  }

  command <<<
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
  >>>

  runtime {
    memory: "~{mem} GB"
    modules: "~{modules}"
    timeout: "~{timeout}"
  }

  output {
    File segments = "~{outputFileNamePrefix}.seg"
    File segmentsWithSubclonalStatus = "~{outputFileNamePrefix}.seg.txt"
    File estimatedCopyNumber = "~{outputFileNamePrefix}.cna.seg"
    File convergedParameters = "~{outputFileNamePrefix}.params.txt"
    File correctedDepth = "~{outputFileNamePrefix}.correctedDepth.txt"
    File rData = "~{outputFileNamePrefix}.RData"
    File plots = "~{outputFileNamePrefix}_plots.tar.gz"
  }

  parameter_meta {
    outputFileNamePrefix: "Output prefix to prefix output file names with."
    wig: "Tumor WIG file."
    normalWig: "Normal WIG file. Default: [NULL]."
    gcWig: "GC-content WIG file."
    mapWig: "Mappability score WIG file. Default: [NULL]."
    normalPanel: "Median corrected depth from panel of normals. Default: [NULL]."
    exonsBed: "Bed file containing exon regions. Default: [NULL]."
    centromere: "File containing Centromere locations; if not provided then will use hg19 version from ichorCNA package."
    minMapScore: "Include bins with a minimum mappability score of this value. Default: [0.9]."
    rmCentromereFlankLength: "Length of region flanking centromere to remove. Default: [1e+05]."
    normal: "Initial normal contamination; can be more than one value if additional normal initializations are desired. Default: [0.5]"
    scStates: "Subclonal states to consider."
    coverage: "PICARD sequencing coverage."
    lambda: "Initial Student's t precision; must contain 4 values (e.g. c(1500,1500,1500,1500)); if not provided then will automatically use based on variance of data."
    lambdaScaleHyperParam: "Hyperparameter (scale) for Gamma prior on Student's-t precision. Default: [3]."
    ploidy: "Initial tumour ploidy; can be more than one value if additional ploidy initializations are desired. Default: [2]"
    maxCN: "Total clonal CN states."
    estimateNormal: "Estimate normal?"
    estimateScPrevalence: "Estimate subclonal prevalence?"
    estimatePloidy: "Estimate tumour ploidy?"
    maxFracCNASubclone: "Exclude solutions with fraction of subclonal events greater than this value. Default: [0.7]."
    maxFracGenomeSubclone: "Exclude solutions with subclonal genome fraction greater than this value. Default: [0.5]."
    minSegmentBins: "Minimum number of bins for largest segment threshold required to estimate tumor fraction; if below this threshold, then will be assigned zero tumor fraction."
    altFracThreshold: "Minimum proportion of bins altered required to estimate tumor fraction; if below this threshold, then will be assigned zero tumor fraction. Default: [0.05]."
    chrNormalize: "Specify chromosomes to normalize GC/mappability biases. Default: [c(1:22)]."
    chrTrain: "Specify chromosomes to estimate params. Default: [c(1:22)]."
    chrs: "Specify chromosomes to analyze."
    genomeBuild: "Genome build. Default: [hg19]."
    genomeStyle: "NCBI or UCSC chromosome naming convention; use UCSC if desired output is to have \"chr\" string. [Default: NCBI]."
    normalizeMaleX: "If male, then normalize chrX by median. Default: [TRUE]."
    fracReadsInChrYForMale: "Threshold for fraction of reads in chrY to assign as male. Default: [0.001]."
    includeHOMD: "If FALSE, then exclude HOMD state. Useful when using large bins (e.g. 1Mb). Default: [FALSE]."
    txnE: "Self-transition probability. Increase to decrease number of segments. Default: [0.9999999]"
    txnStrength: "Transition pseudo-counts. Exponent should be the same as the number of decimal places of --txnE. Default: [1e+07]."
    plotFileType: "File format for output plots. Default: [pdf]."
    plotYLim: "ylim to use for chromosome plots. Default: [c(-2,2)]."
    outDir: "Output Directory. Default: [./]."
    libdir: "Script library path."
    modules: "Environment module name and version to load (space separated) before command execution."
    mem: "Memory (in GB) to allocate to the job."
    timeout: "Maximum amount of time (in hours) the task can run for."
  }

  meta {
    output_meta: {
      segments: "Segments called by the Viterbi algorithm.  Format is compatible with IGV.",
      segmentsWithSubclonalStatus: "Same as `segments` but also includes subclonal status of segments (0=clonal, 1=subclonal). Format not compatible with IGV.",
      estimatedCopyNumber: "Estimated copy number, log ratio, and subclone status for each bin/window.",
      convergedParameters: "Final converged parameters for optimal solution. Also contains table of converged parameters for all solutions.",
      correctedDepth: "Log2 ratio of each bin/window after correction for GC and mappability biases.",
      rData: "Saved R image after ichorCNA has finished. Results for all solutions will be included.",
      plots: "Archived directory of plots."
    }
  }
}
