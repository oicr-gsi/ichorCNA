version 1.0

import "imports/pull_bwa.wdl" as bwaMem

struct InputGroup {
  File fastqR1
  File fastqR2
  String readGroups
}

struct PdfOutput {
  Array[Pair[File,Map[String,String]]]+ pdfs
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

    call preMergeBamMetrics as preMergeBamMetricsFastqInput{
      input:
        bam = bwaMem.bwaMemBam,
        outputFileNamePrefix = outputFileNamePrefix
    }

    if (length(inputGroups_) > 1 ) {
      call bamMerge {
        input:
          bams = bwaMem.bwaMemBam,
          outputFileNamePrefix = outputFileNamePrefix
      }
      #get reads of each input file
    }

    if (length(inputGroups_) == 1 ) {
      File bwaMemBam = bwaMem.bwaMemBam[0]
    }

  }

  if(inputType=="bam" && defined(inputBam)){
    Array[File] inputBam_ = select_first([inputBam,[]])

    call preMergeBamMetrics {
      input:
        bam = inputBam_,
        outputFileNamePrefix = outputFileNamePrefix
    }

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

  call getMetrics {
    input:
      inputbam = select_first([bamMerge.outputMergedBam,bwaMemBam,inputBamMerge.outputMergedBam,singleInputBam]),
      params = runIchorCNA.convergedParameters,
      outputFileNamePrefix = outputFileNamePrefix
  }

  call createJson {
    input:
      bamMetrics = getMetrics.bamMetrics,
      preBamMetrics = select_first([preMergeBamMetrics.preMergeMetrics,preMergeBamMetricsFastqInput.preMergeMetrics]),
      allSolsMetrics = getMetrics.all_sols_metrics,
      plotsFile = runIchorCNA.plotsTxt,
      outputFileNamePrefix = outputFileNamePrefix
  }

  output {
    Pair[File,Map[String,String]] solution1 = createJson.pdfOutput.pdfs[0]
    Pair[File,Map[String,String]] solution2 = createJson.pdfOutput.pdfs[1]
    Pair[File,Map[String,String]] solution3 = createJson.pdfOutput.pdfs[2]
    Pair[File,Map[String,String]] solution4 = createJson.pdfOutput.pdfs[3]
    Pair[File,Map[String,String]] solution5 = createJson.pdfOutput.pdfs[4]
    Pair[File,Map[String,String]] solution6 = createJson.pdfOutput.pdfs[5]
    Pair[File,Map[String,String]] solution7 = createJson.pdfOutput.pdfs[6]
    Pair[File,Map[String,String]] solution8 = createJson.pdfOutput.pdfs[7]
    Pair[File,Map[String,String]] solution9 = createJson.pdfOutput.pdfs[8]
    Pair[File,Map[String,String]] solution10 = createJson.pdfOutput.pdfs[9]
    Pair[File,Map[String,String]] solution11 = createJson.pdfOutput.pdfs[10]
    Pair[File,Map[String,String]] solution12 = createJson.pdfOutput.pdfs[11]
    Pair[File,Map[String,String]] solution13 = createJson.pdfOutput.pdfs[12]
    Pair[File,Map[String,String]] solution14 = createJson.pdfOutput.pdfs[13]
    Pair[File,Map[String,String]] solution15 = createJson.pdfOutput.pdfs[14]
    Pair[File,Map[String,String]] solution16 = createJson.pdfOutput.pdfs[15]
    File? bam = indexBam.outbam
    File? bamIndex = indexBam.bamIndex
    File jsonMetrics = createJson.metricsJson
    File segments = runIchorCNA.segments
    File segmentsWithSubclonalStatus = runIchorCNA.segmentsWithSubclonalStatus
    File estimatedCopyNumber = runIchorCNA.estimatedCopyNumber
    File convergedParameters = runIchorCNA.convergedParameters
    File correctedDepth = runIchorCNA.correctedDepth
    File rData = runIchorCNA.rData
    File plots = runIchorCNA.plots
  }

  meta {
    author: "Beatriz Lujan Toro"
    email: "beatriz.lujantoro@oicr.on.ca"
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
      },
      {
        name: "python/3.9",
        url: "python.org"
      },
      {
        name: "pandas/1.4.2",
        url: "https://pandas.pydata.org/"
      }
    ]
    output_meta: {
      pdf: "Annotations for pdf files produced by ichorCNA, each pdf produced is annotated with tumor fraction, ploidy and log likelihood.",
      bam: "Bam file used as input to ichorCNA (only produced when provisionBam is True)",
      bamIndex: "Bam index for bam file used as input to ichorCNA (only produced when provisionBam is True)",
      jsonMetrics: "Report on bam coverage, read counts and ichorCNA metrics.",
      segments: "Segments called by the Viterbi algorithm.  Format is compatible with IGV.",
      segmentsWithSubclonalStatus: "Same as segments but also includes subclonal status of segments (0=clonal, 1=subclonal). Format not compatible with IGV.",
      estimatedCopyNumber: "Estimated copy number, log ratio, and subclone status for each bin/window.",
      convergedParameters: "Final converged parameters for optimal solution. Also contains table of converged parameters for all solutions.",
      correctedDepth: "Log2 ratio of each bin/window after correction for GC and mappability biases.",
      rData: "Saved R image after ichorCNA has finished. Results for all solutions will be included.",
      plots: "Archived directory of plots.",
      solution1: "Plots for solution 1.",
      solution2: "Plots for solution 2.",
      solution3: "Plots for solution 3.",
      solution4: "Plots for solution 4.",
      solution5: "Plots for solution 5.",
      solution6: "Plots for solution 6.",
      solution7: "Plots for solution 7.",
      solution8: "Plots for solution 8.",
      solution9: "Plots for solution 9.",
      solution10: "Plots for solution 10.",
      solution11: "Plots for solution 11.",
      solution12: "Plots for solution 12.",
      solution13: "Plots for solution 13.",
      solution14: "Plots for solution 14.",
      solution15: "Plots for solution 15.",
      solution16: "Plots for solution 16."

    }
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

task preMergeBamMetrics {
  input {
    Array[File] bam
    String outputFileNamePrefix
    Int jobMemory = 8
    String modules = "samtools/1.14"
    Int timeout = 12
  }

  parameter_meta {
    bam: "Input bam pre merge."
    outputFileNamePrefix: "Output prefix to prefix output file names with."
    jobMemory: "Memory (in GB) to allocate to the job."
    modules: "Environment module name and version to load (space separated) before command execution."
    timeout: "Maximum amount of time (in hours) the task can run for."
  }

  command <<<
  set -exo pipefail

  echo run_lane,read_count > ~{outputFileNamePrefix}_pre_merge_bam_metrics.csv
  for file in ~{sep=' ' bam}
  do
    run=$(samtools view -H "${file}" | grep '^@RG' | cut -f 2 | cut -f 2 -d ":")
    read_count=$(samtools view -c "${file}")
    echo $run,$read_count >> ~{outputFileNamePrefix}_pre_merge_bam_metrics.csv
  done;

  >>>

  output {
  File preMergeMetrics = "~{outputFileNamePrefix}_pre_merge_bam_metrics.csv"
  }

  runtime {
    memory: "~{jobMemory} GB"
    modules: "~{modules}"
    timeout: "~{timeout}"
  }

  meta {
    output_meta: {
      preMergeMetrics: "csv file with bam metrics."
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

  parameter_meta {
    inputbam: "Input bam."
    jobMemory: "Memory (in GB) to allocate to the job."
    modules: "Environment module name and version to load (space separated) before command execution."
    timeout: "Maximum amount of time (in hours) the task can run for."
  }

  String resultBai = "~{basename(inputbam)}.bai"

  command <<<
  samtools index ~{inputbam} ~{resultBai}
  >>>

  runtime {
    memory: "~{jobMemory} GB"
    modules: "~{modules}"
    timeout: "~{timeout}"
  }

  output {
    File outbam = "~{inputbam}"
    File bamIndex = "~{resultBai}"
  }

  meta {
    output_meta: {
      outbam: "alignment file in bam format used for the analysis (merged if input is multiple fastqs or bams).",
      bamIndex: "output index file for bam aligned to genome."
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

    #create txt file with plot full path
    ls $PWD/~{outputFileNamePrefix}/*genomeWide_n* > "~{outputFileNamePrefix}"_plots.txt
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
    File plotsTxt = "~{outputFileNamePrefix}_plots.txt"
    File solution1 = "~{outputFileNamePrefix}/~{outputFileNamePrefix}_genomeWide_n0.2-p2.pdf"
    File solution2 = "~{outputFileNamePrefix}/~{outputFileNamePrefix}_genomeWide_n0.2-p3.pdf"
    File solution3 = "~{outputFileNamePrefix}/~{outputFileNamePrefix}_genomeWide_n0.3-p2.pdf"
    File solution4 = "~{outputFileNamePrefix}/~{outputFileNamePrefix}_genomeWide_n0.3-p3.pdf"
    File solution5 = "~{outputFileNamePrefix}/~{outputFileNamePrefix}_genomeWide_n0.4-p2.pdf"
    File solution6 = "~{outputFileNamePrefix}/~{outputFileNamePrefix}_genomeWide_n0.4-p3.pdf"
    File solution7 = "~{outputFileNamePrefix}/~{outputFileNamePrefix}_genomeWide_n0.5-p2.pdf"
    File solution8 = "~{outputFileNamePrefix}/~{outputFileNamePrefix}_genomeWide_n0.5-p3.pdf"
    File solution9 = "~{outputFileNamePrefix}/~{outputFileNamePrefix}_genomeWide_n0.6-p2.pdf"
    File solution10 = "~{outputFileNamePrefix}/~{outputFileNamePrefix}_genomeWide_n0.6-p3.pdf"
    File solution11 = "~{outputFileNamePrefix}/~{outputFileNamePrefix}_genomeWide_n0.7-p2.pdf"
    File solution12 = "~{outputFileNamePrefix}/~{outputFileNamePrefix}_genomeWide_n0.7-p3.pdf"
    File solution13 = "~{outputFileNamePrefix}/~{outputFileNamePrefix}_genomeWide_n0.8-p2.pdf"
    File solution14 = "~{outputFileNamePrefix}/~{outputFileNamePrefix}_genomeWide_n0.8-p3.pdf"
    File solution15 = "~{outputFileNamePrefix}/~{outputFileNamePrefix}_genomeWide_n0.9-p2.pdf"
    File solution16 = "~{outputFileNamePrefix}/~{outputFileNamePrefix}_genomeWide_n0.9-p3.pdf"

  }


  meta {
    output_meta: {
      segments: "Segments called by the Viterbi algorithm.  Format is compatible with IGV.",
      segmentsWithSubclonalStatus: "Same as `segments` but also includes subclonal status of segments (0=clonal, 1=subclonal). Format not compatible with IGV.",
      estimatedCopyNumber: "Estimated copy number, log ratio, and subclone status for each bin/window.",
      convergedParameters: "Final converged parameters for optimal solution. Also contains table of converged parameters for all solutions.",
      correctedDepth: "Log2 ratio of each bin/window after correction for GC and mappability biases.",
      rData: "Saved R image after ichorCNA has finished. Results for all solutions will be included.",
      plots: "Archived directory of plots.",
      plotsTxt: "Text file with the full path to the solution 1-16 pdfs.",
      solution1: "Plots for solution 1.",
      solution2: "Plots for solution 2.",
      solution3: "Plots for solution 3.",
      solution4: "Plots for solution 4.",
      solution5: "Plots for solution 5.",
      solution6: "Plots for solution 6.",
      solution7: "Plots for solution 7.",
      solution8: "Plots for solution 8.",
      solution9: "Plots for solution 9.",
      solution10: "Plots for solution 10.",
      solution11: "Plots for solution 11.",
      solution12: "Plots for solution 12.",
      solution13: "Plots for solution 13.",
      solution14: "Plots for solution 14.",
      solution15: "Plots for solution 15.",
      solution16: "Plots for solution 16."
    }
  }
}

task getMetrics {
  input {
    File inputbam
    File params
    String outputFileNamePrefix
    Int jobMemory = 8
    String modules = "samtools/1.14"
    Int timeout = 12
  }

  parameter_meta {
    inputbam: "Input bam."
    outputFileNamePrefix: "Output prefix to prefix output file names with."
    jobMemory: "Memory (in GB) to allocate to the job."
    modules: "Environment module name and version to load (space separated) before command execution."
    timeout: "Maximum amount of time (in hours) the task can run for."
  }

  command <<<
  echo coverage,read_count,tumor_fraction,ploidy > ~{outputFileNamePrefix}_bam_metrics.csv
  coverage=$(samtools coverage ~{inputbam} | grep -P "^chr\d+\t|^chrX\t|^chrY\t" | awk '{ space += ($3-$2)+1; bases += $7*($3-$2);} END { print bases/space }')
  read_count=$(samtools stats ~{inputbam} | head -n 8 | tail -n 1 | cut -f 3)
  tumor_fraction=$(cat ~{params} | head -n 2 | tail -n 1 | cut -f 2)
  ploidy=$(cat ~{params} | head -n 2 | tail -n 1 | cut -f 3)
  echo $coverage,$read_count,$tumor_fraction,$ploidy >> ~{outputFileNamePrefix}_bam_metrics.csv
  cat ~{params} | tail -n 17 > ~{outputFileNamePrefix}_all_sols_metrics.csv
  >>>

  output {
  File bamMetrics = "~{outputFileNamePrefix}_bam_metrics.csv"
  File all_sols_metrics = "~{outputFileNamePrefix}_all_sols_metrics.csv"
  }

  runtime {
    memory: "~{jobMemory} GB"
    modules: "~{modules}"
    timeout: "~{timeout}"
  }

  meta {
    output_meta: {
      bamMetrics: "Metrics collected from bam file used for ichorCNA, to be used as input for final json metrics collection (createJson task).",
      all_sols_metrics: "Collected metrics from each solution stored in the params file, to be used as input for final json metrics collection (createJson task)."
    }
  }
}

task createJson {
  input {
    File preBamMetrics
    File bamMetrics
    File allSolsMetrics
    File plotsFile
    String outputFileNamePrefix
    Int jobMemory = 8
    String modules = "pandas/1.4.2"
    Int timeout = 12
  }

  parameter_meta {
    preBamMetrics: "pre-merge bam metrics."
    bamMetrics: "bam metrics."
    allSolsMetrics: "metrics for all solutions"
    outputFileNamePrefix: "Output prefix to prefix output file names with."
    jobMemory: "Memory (in GB) to allocate to the job."
    modules: "Environment module name and version to load (space separated) before command execution."
    timeout: "Maximum amount of time (in hours) the task can run for."
  }

  command <<<
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
      len_pdf_sol = len(line.split("_")[-1])
      pdf_solution = line.split("_")[-1][0:(len_pdf_sol-4)]
      pdf_dict["left"] = line
      pdf_dict["right"] = {}
      pdf_dict["right"] = metrics_dict["solutions"][pdf_solution]
      output_list.append(pdf_dict)
    output_dict = {}
    output_dict["pdfs"] = output_list

    with open("~{outputFileNamePrefix}_outputs.json", "w") as outPdfJson:
        json.dump(output_dict, outPdfJson)

    CODE
  >>>

  output {
    File metricsJson = "~{outputFileNamePrefix}_metrics.json"
    PdfOutput pdfOutput = read_json("~{outputFileNamePrefix}_outputs.json")
  }

  runtime {
    memory: "~{jobMemory} GB"
    modules: "~{modules}"
    timeout: "~{timeout}"
  }

  meta {
    output_meta: {
      metricsJson: "json file reporting mean coverage, total reads, lanes sequenced, reads per lane as well as ichorCNA reported ploidy, tumor_fraction for selected and all reported solutions.",
      out: "Annotated output."
    }
  }
}
