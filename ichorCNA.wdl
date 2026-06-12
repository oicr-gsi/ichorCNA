version 1.0

import "imports/pull_bamQC.wdl" as bamQC

struct ichorCNAResources {
    String gcWig
    String mapWig
    String normalPanel
    String centromere
    String refFasta
}

struct PdfOutput {
  Array[Pair[File,Map[String,String]]]+ pdfs
}


workflow ichorCNA {
  input {
    Array[File]? inputBam
    File? inputCram
    String outputFileNamePrefix
    Int windowSize
    Int minimumMappingQuality
    String chromosomesToAnalyze
    Boolean provisionBam
    String reference
  }

  parameter_meta {
    inputBam: "Array of one or multiple bam files. Provide either inputBam or inputCram (not both). BAM input uses the readCounter wig-generation path."
    inputCram: "Single cram file. Provide either inputBam or inputCram (not both). CRAM input uses the mosdepth wig-generation path."
    outputFileNamePrefix: "Output prefix to prefix output file names with."
    windowSize: "The size of non-overlapping windows."
    minimumMappingQuality: "Mapping quality value below which reads are ignored."
    chromosomesToAnalyze: "Chromosomes in the bam reference file."
    provisionBam: "Boolean, to provision out bam file and coverage metrics"
    reference: "The genome reference build. for example: hg19, hg38"
  }


  Map[String,ichorCNAResources] resources = {
    "hg19": {
      "gcWig": "$ICHORCNA_ROOT/lib/R/library/ichorCNA/extdata/gc_hg19_1000kb.wig",
      "mapWig": "$ICHORCNA_ROOT/lib/R/library/ichorCNA/extdata/map_hg19_1000kb.wig",
      "normalPanel": "$ICHORCNA_ROOT/lib/R/library/ichorCNA/extdata/HD_ULP_PoN_1Mb_median_normAutosome_mapScoreFiltered_median.rds",
      "centromere": "$ICHORCNA_ROOT/lib/R/library/ichorCNA/extdata/GRCh37.p13_centromere_UCSC-gapTable.txt",
      "refFasta": "/.mounts/labs/gsi/modulator/sw/data/hg19-p13/hg19_random.fa"
    },
    "hg38": {
      "gcWig": "$ICHORCNA_ROOT/lib/R/library/ichorCNA/extdata/gc_hg38_1000kb.wig",
      "mapWig": "$ICHORCNA_ROOT/lib/R/library/ichorCNA/extdata/map_hg38_1000kb.wig",
      "normalPanel": "$ICHORCNA_ROOT/lib/R/library/ichorCNA/extdata/HD_ULP_PoN_hg38_1Mb_median_normAutosome_median.rds",
      "centromere": "$ICHORCNA_ROOT/lib/R/library/ichorCNA/extdata/GRCh38.GCA_000001405.2_centromere_acen.txt",
      "refFasta": "$HG38_ULTIMA_ROOT/Homo_sapiens_assembly38.fasta"
    },
    "hg38_noAlt": {
      "gcWig": "$ICHORCNA_ROOT/lib/R/library/ichorCNA/extdata/gc_hg38_1000kb.wig",
      "mapWig": "$ICHORCNA_ROOT/lib/R/library/ichorCNA/extdata/map_hg38_1000kb.wig",
      "normalPanel": "$ICHORCNA_ROOT/lib/R/library/ichorCNA/extdata/HD_ULP_PoN_hg38_1Mb_median_normAutosome_median.rds",
      "centromere": "$ICHORCNA_ROOT/lib/R/library/ichorCNA/extdata/GRCh38.GCA_000001405.2_centromere_acen.txt",
      "refFasta": "/.mounts/labs/gsi/modulator/sw/data/hg38-p12/hg38_random.fa"
    }
  }

  # BAM input -> readCounter wig-generation path
  if (defined(inputBam)) {
    Array[File] inputBam_ = select_first([inputBam])

    call preMergeBamMetrics as preMergeBamMetricsBam {
      input:
        bam = inputBam_,
        outputFileNamePrefix = outputFileNamePrefix
    }

    if (length(inputBam_) > 1 ) {
      call bamMerge {
        input:
          bams = inputBam_,
          outputFileNamePrefix = outputFileNamePrefix
      }
    }
    File finalBam = select_first([
      bamMerge.outputMergedBam,
      inputBam_[0]
    ])

    if(provisionBam){
      call indexBam {
        input:
          inputbam = finalBam
      }
    }

    call runReadCounter{
      input:
        bam= finalBam,
        outputFileNamePrefix=outputFileNamePrefix,
        windowSize=windowSize,
        minimumMappingQuality=minimumMappingQuality,
        chromosomesToAnalyze=chromosomesToAnalyze
    }

    call bamQC.bamQC {
      input:
          bamFile = finalBam,
          outputFileNamePrefix = outputFileNamePrefix
    }
  }

  # CRAM input -> mosdepth wig-generation path
  if (defined(inputCram)) {
    File inputCram_ = select_first([inputCram])
    String cramRefFasta = resources [ reference ].refFasta

    call runMosdepth {
      input:
        cram = inputCram_,
        refFasta = cramRefFasta,
        outputFileNamePrefix = outputFileNamePrefix,
        windowSize = windowSize
    }
  }

  File wig = select_first([runReadCounter.wig, runMosdepth.wig])
  Array[String] chrs = select_first([runReadCounter.ichorCNAchrs, runMosdepth.ichorCNAchrs])

  call runIchorCNA {
    input:
      outputFileNamePrefix=outputFileNamePrefix,
      chrs=chrs,
      wig=wig,
      gcWig = resources [ reference ].gcWig,
      mapWig = resources [ reference ].mapWig,
      normalPanel = resources [ reference ].normalPanel,
      centromere = resources [ reference ].centromere,
      genomeBuild=reference
  }

  # BAM path: per-lane counts come from preMergeBamMetricsBam (pre-merge), and
  # getMetrics computes coverage/reads on the merged bam.
  if (defined(inputBam)) {
    call getMetrics {
      input:
        inputbam = select_first([finalBam]),
        params = runIchorCNA.convergedParameters,
        outputFileNamePrefix = outputFileNamePrefix
    }
  }

  # CRAM path: a single cram is one lane, so one task does the (expensive) cram
  # decode once and emits both the lane-level and sample-level metrics.
  if (defined(inputCram)) {
    call getCramMetrics {
      input:
        cram = select_first([inputCram]),
        refFasta = select_first([cramRefFasta]),
        params = runIchorCNA.convergedParameters,
        outputFileNamePrefix = outputFileNamePrefix
    }
  }

  File preBamMetrics = select_first([preMergeBamMetricsBam.preMergeMetrics, getCramMetrics.laneMetrics])
  File bamMetrics = select_first([getMetrics.bamMetrics, getCramMetrics.bamMetrics])
  File allSolsMetrics = select_first([getMetrics.all_sols_metrics, getCramMetrics.all_sols_metrics])

  call createJson {
    input:
      bamMetrics = bamMetrics,
      preBamMetrics = preBamMetrics,
      allSolsMetrics = allSolsMetrics,
      plotsFile = runIchorCNA.plotsTxt,
      outputFileNamePrefix = outputFileNamePrefix
  }

  output {
    Pair[File,Map[String,String]] genomeWideAll = createJson.pdfOutput.pdfs[0]
    Pair[File,Map[String,String]] genomeWide = createJson.pdfOutput.pdfs[1]
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
    File? bamQCresult = bamQC.result
  }

  meta {
    author: "Beatriz Lujan Toro and Aditi Nallan"
    email: "beatriz.lujantoro@oicr.on.ca and anallan@oicr.on.ca"
    description: "Workflow for estimating the fraction of tumor in cell-free DNA from sWGS (shallow Whole Genome Sequencing). ichorCNA can be used to inform the presence or absence of tumor-derived DNA and to guide the decision to perform whole exome or deeper whole genome sequencing. Furthermore, the quantitative estimate of tumor fraction can we used to calibrate the desired depth of sequencing to reach statistical power for identifying mutations in cell-free DNA. Finally, ichorCNA can be use to detect large-scale copy number alterations from large cohorts by taking advantage of the cost-effective approach of ultra-low-pass sequencing."
    dependencies: [
      {
        name: "samtools/1.14",
        url: "http://www.htslib.org/"
      },
      {
        name: "hmmcopy-utils/0.1.1",
        url: "https://shahlab.ca/projects/hmmcopy_utils/"
      },
      {
        name: "mosdepth/0.3.8",
        url: "https://github.com/brentp/mosdepth"
      },
      {
        name: "ichorcna/0.2",
        url: "https://github.com/broadinstitute/ichorCNA"
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
    pdf: {
        description: "Annotations for pdf files produced by ichorCNA, each pdf is annotated with the tumor fraction, ploidy for the selected solution.",
        vidarr_label: "pdf"
    },
    bam: {
        description: "Bam file used as input to ichorCNA",
        vidarr_label: "bam"
    },
    bamIndex: {
        description: "Bam index for bam file used as input to ichorCNA",
        vidarr_label: "bamIndex"
    },
    bamQCresult: {
        description: "bamQC report.",
        vidarr_label: "bamQCresult"
    },
    jsonMetrics: {
        description: "Report on bam coverage, read counts and ichorCNA metrics.",
        vidarr_label: "jsonMetrics"
    },
    segments: {
        description: "Segments called by the Viterbi algorithm.  Format is compatible with IGV.",
        vidarr_label: "segments"
    },
    segmentsWithSubclonalStatus: {
        description: "Same as segments but also includes subclonal status of segments (0=clonal, 1=subclonal). Format not compatible with IGV.",
        vidarr_label: "segmentsWithSubclonalStatus"
    },
    estimatedCopyNumber: {
        description: "Estimated copy number, log ratio, and subclone status for each bin/window.",
        vidarr_label: "estimatedCopyNumber"
    },
    convergedParameters: {
        description: "Final converged parameters for optimal solution. Also contains table of converged parameters for all solutions.",
        vidarr_label: "convergedParameters"
    },
    correctedDepth: {
        description: "Log2 ratio of each bin/window after correction for GC and mappability biases.",
        vidarr_label: "correctedDepth"
    },
    rData: {
        description: "Saved R image after ichorCNA has finished. Results for all solutions will be included.",
        vidarr_label: "rData"
    },
    plots: {
        description: "Archived directory of plots.",
        vidarr_label: "plots"
    },
    genomeWideAll: "Genome wide plots for each solution",
    genomeWide: "Genome wide plots for the selected solution"
  }
  }
}

task bamMerge{
  input {
    Array[File] bams
    String outputFileNamePrefix
    Int   jobMemory = 32
    String modules  = "samtools/1.14"
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
    String? refFasta
    String outputFileNamePrefix
    Int jobMemory = 8
    String modules = "samtools/1.14"
    Int timeout = 12
  }

  parameter_meta {
    bam: "Input bam (or cram) pre merge."
    refFasta: "Reference FASTA, required when the input is cram so samtools can decode it. Default: [NULL]."
    outputFileNamePrefix: "Output prefix to prefix output file names with."
    jobMemory: "Memory (in GB) to allocate to the job."
    modules: "Environment module name and version to load (space separated) before command execution."
    timeout: "Maximum amount of time (in hours) the task can run for."
  }

  command <<<
  set -euo pipefail

  echo run,read_count > ~{outputFileNamePrefix}_pre_merge_bam_metrics.csv
  for file in ~{sep=' ' bam}
  do
    run=$(samtools view ~{"--reference " + refFasta} -H "${file}" | grep '^@RG' | cut -f 2 | cut -f 2 -d ":" | cut -f 1 -d "-")
    run="${run%%$'\n'*}"   # keep only the first run name if the file has multiple @RG lines
    read_count=$(samtools stats ~{"--reference " + refFasta} "${file}" | grep ^SN | grep "raw total sequences" | cut -f 3)
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
    Int jobMemory = 12
    String modules = "samtools/1.14"
    Int timeout = 48
  }

  parameter_meta {
    inputbam: "Input bam."
    jobMemory: "Memory (in GB) to allocate to the job."
    modules: "Environment module name and version to load (space separated) before command execution."
    timeout: "Maximum amount of time (in hours) the task can run for."
  }

  String resultBai = "~{basename(inputbam)}.bai"

  command <<<
  set -euo pipefail
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
      outbam: "alignment file in bam format used for the analysis (merged if input is multiple bams).",
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
    String modules = "samtools/1.14 hmmcopy-utils/0.1.1"
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
    set -euo pipefail

    samtools index ~{bam}

    # calculate chromosomes to analyze (with reads) from input data
    CHROMOSOMES_WITH_READS=$(samtools idxstats ~{bam} | awk '$3 > 0' - | cut -f1 | grep -Ew $(tr ',' '|' <<<  '~{chromosomesToAnalyze}') | paste -s -d, -)

    # write out a chromosomes with reads for ichorCNA
    # split onto new lines (for wdl read_lines), exclude chrY, remove chr prefix, wrap in single quotes for ichorCNA
    echo "${CHROMOSOMES_WITH_READS}" | tr ',' '\n' | grep -v chrY | sed "s/chr//g" | sort -V | sed -e "s/\(.*\)/'\1'/" > ichorCNAchrs.txt

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

task runMosdepth {
  input{
    File cram
    String refFasta
    String outputFileNamePrefix
    Int windowSize
    Int threads = 8
    Int mem = 16
    String modules = "mosdepth/0.3.8 samtools/1.14 hg38-ultima/v0"
    Int timeout = 12
  }

  parameter_meta {
    cram: "Input cram."
    refFasta: "Reference FASTA used to decode the cram and by mosdepth (--fasta)."
    outputFileNamePrefix: "Output prefix to prefix output file names with."
    windowSize: "The size of non-overlapping windows (mosdepth --by)."
    threads: "Number of threads for mosdepth."
    mem: "Memory (in GB) to allocate to the job."
    modules: "Environment module name and version to load (space separated) before command execution."
    timeout: "Maximum amount of time (in hours) the task can run for."
  }

  command <<<
    set -euo pipefail

    mkdir -p mosdepth

    samtools index ~{cram}

    # per-window mean coverage with mosdepth
    mosdepth \
    -t ~{threads} \
    --fast-mode \
    --by ~{windowSize} \
    --fasta ~{refFasta} \
    mosdepth/~{outputFileNamePrefix} \
    ~{cram}

    # convert the binned coverage to a fixedStep WIG; strip the "chr" prefix from
    # chrom names so they match the gcWig/mapWig reference WIGs (as the bam path does)
    zcat mosdepth/~{outputFileNamePrefix}.regions.bed.gz | \
    awk 'BEGIN { OFS="\t" }
    $1 ~ /^chr([1-9]|1[0-9]|2[0-2]|X|Y)$/ {
      chr = $1; start = $2; end = $3; cov = $4;
      if (NR == 1 || chr != prev_chr || start != prev_end) {
        print "fixedStep chrom=" chr " start=" (start+1) " step=" (end - start) " span=" (end - start);
      }
      print cov;
      prev_chr = chr; prev_end = end;
    }' | sed "s/chrom=chr/chrom=/" > ~{outputFileNamePrefix}.wig

    # write out chromosomes with data for ichorCNA
    # (chr prefix already stripped above, exclude Y, sort, wrap in single quotes)
    grep -E "^fixedStep" ~{outputFileNamePrefix}.wig | \
    sed -E 's/.*chrom=([^ ]+).*/\1/' | \
    grep -vw Y | sort -uV | sed -e "s/\(.*\)/'\1'/" > ichorCNAchrs.txt
  >>>

  runtime {
    memory: "~{mem} GB"
    cpu: "~{threads}"
    modules: "~{modules}"
    timeout: "~{timeout}"
  }

  output {
    File wig = "~{outputFileNamePrefix}.wig"
    Array[String] ichorCNAchrs = read_lines("ichorCNAchrs.txt")
  }

  meta {
    output_meta: {
      wig: "Per-window coverage file in WIG format",
      ichorCNAchrs: "Chromosomes with data for ichorCNA (\"chr\" stripped from the name)"
    }
  }
}

task runIchorCNA {
  input {
    String outputFileNamePrefix
    File wig
    File? normalWig
    String gcWig
    String mapWig
    String normalPanel
    String? exonsBed
    String centromere
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
    String genomeBuild
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
    set -euo pipefail

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
    ls $PWD/~{outputFileNamePrefix}/*genomeWide_all_sols.pdf > "~{outputFileNamePrefix}"_plots.txt
    ls $PWD/~{outputFileNamePrefix}/*genomeWide.pdf >> "~{outputFileNamePrefix}"_plots.txt
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
    File genomeWideAll = "~{outputFileNamePrefix}/~{outputFileNamePrefix}_genomeWide_all_sols.pdf"
    File genomeWide = "~{outputFileNamePrefix}/~{outputFileNamePrefix}_genomeWide.pdf"
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
      genomeWideAll: "Genome wide plots for each solution",
      genomeWide: "Genome wide plots for the selected solution"
    }
  }
}

task getMetrics {
  input {
    File inputbam
    String? refFasta
    File params
    String outputFileNamePrefix
    Int jobMemory = 8
    String modules = "samtools/1.14"
    Int timeout = 12
  }

  parameter_meta {
    inputbam: "Input bam (or cram)."
    refFasta: "Reference FASTA, required when the input is cram so samtools can decode it. Default: [NULL]."
    outputFileNamePrefix: "Output prefix to prefix output file names with."
    jobMemory: "Memory (in GB) to allocate to the job."
    modules: "Environment module name and version to load (space separated) before command execution."
    timeout: "Maximum amount of time (in hours) the task can run for."
  }

  command <<<
  set -euo pipefail

  echo coverage,read_count,tumor_fraction,ploidy > ~{outputFileNamePrefix}_bam_metrics.csv
  coverage=$(samtools coverage ~{"--reference " + refFasta} ~{inputbam} | grep -P "^chr\d+\t|^chrX\t|^chrY\t" | awk '{ space += ($3-$2)+1; bases += $7*($3-$2);} END { print bases/space }')
  read_count=$(samtools stats ~{"--reference " + refFasta} ~{inputbam} | grep ^SN | grep "raw total sequences" | cut -f 3)
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

task getCramMetrics {
  input {
    File cram
    String refFasta
    File params
    String outputFileNamePrefix
    Int jobMemory = 8
    String modules = "samtools/1.14 hg38-ultima/v0"
    Int timeout = 24
  }

  parameter_meta {
    cram: "Input cram."
    refFasta: "Reference FASTA used by samtools to decode the cram."
    params: "Converged parameters file from runIchorCNA (tumor fraction / ploidy)."
    outputFileNamePrefix: "Output prefix to prefix output file names with."
    jobMemory: "Memory (in GB) to allocate to the job."
    modules: "Environment module name and version to load (space separated) before command execution."
    timeout: "Maximum amount of time (in hours) the task can run for."
  }

  # Decoding a cram is expensive, so this runs once and emits both the
  # lane-level (preBamMetrics) and sample-level (bamMetrics) CSVs, replacing the
  # separate preMergeBamMetrics + getMetrics passes the bam path uses. A single
  # cram is one "lane" (nothing is merged).
  command <<<
  set -euo pipefail

  run=$(samtools view --reference ~{refFasta} -H ~{cram} | grep '^@RG' | cut -f 2 | cut -f 2 -d ":" | cut -f 1 -d "-")
  run="${run%%$'\n'*}"   # keep only the first run name if the cram has multiple @RG lines
  read_count=$(samtools stats --reference ~{refFasta} ~{cram} | grep ^SN | grep "raw total sequences" | cut -f 3)
  coverage=$(samtools coverage --reference ~{refFasta} ~{cram} | grep -P "^chr\d+\t|^chrX\t|^chrY\t" | awk '{ space += ($3-$2)+1; bases += $7*($3-$2);} END { print bases/space }')
  tumor_fraction=$(cat ~{params} | head -n 2 | tail -n 1 | cut -f 2)
  ploidy=$(cat ~{params} | head -n 2 | tail -n 1 | cut -f 3)

  # lane-level CSV (one row for the single cram) — schema matches preMergeBamMetrics
  echo run,read_count > ~{outputFileNamePrefix}_pre_merge_bam_metrics.csv
  echo $run,$read_count >> ~{outputFileNamePrefix}_pre_merge_bam_metrics.csv

  # sample-level CSV — schema matches getMetrics
  echo coverage,read_count,tumor_fraction,ploidy > ~{outputFileNamePrefix}_bam_metrics.csv
  echo $coverage,$read_count,$tumor_fraction,$ploidy >> ~{outputFileNamePrefix}_bam_metrics.csv

  cat ~{params} | tail -n 17 > ~{outputFileNamePrefix}_all_sols_metrics.csv
  >>>

  output {
  File laneMetrics = "~{outputFileNamePrefix}_pre_merge_bam_metrics.csv"
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
      laneMetrics: "Per-lane read counts (single row for a single cram); fed to createJson as preBamMetrics for lanes_sequenced / reads_per_lane.",
      bamMetrics: "Sample-level coverage and read count plus ichorCNA tumor fraction / ploidy, for final json metrics collection (createJson task).",
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
    set -euo pipefail

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
