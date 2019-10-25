version 1.0

workflow ichorCNA {

  input {
    File bam
    File bamIndex
    String? outputFileNamePrefix
    Int windowSize
    Int minimumMappingQuality
    String chromosomesToAnalyze
  }

  call runReadCounter{
    input:
      bam=bam,
      bamIndex=bamIndex,
      outputFileNamePrefix=outputFileNamePrefix,
      windowSize=windowSize,
      minimumMappingQuality=minimumMappingQuality,
      chromosomesToAnalyze=chromosomesToAnalyze
 }

  call runIchorCNA {
    input:
      outputFileNamePrefix=outputFileNamePrefix,
      wig=runReadCounter.wig
  }

  output {
    File segments = runIchorCNA.segments
    File segmentsWithSubclonalStatus = runIchorCNA.segmentsWithSubclonalStatus
    File estimatedCopyNumber = runIchorCNA.estimatedCopyNumber
    File convergedParameters = runIchorCNA.convergedParameters
    File correctedDepth = runIchorCNA.correctedDepth
    File rData = runIchorCNA.rData
    File plots = runIchorCNA.plots
  }
}

# First step : read counter
task runReadCounter {
  input{
    File bam
    File bamIndex
    String? outputFileNamePrefix = "output"
    Int windowSize
    Int minimumMappingQuality
    String chromosomesToAnalyze
    Int? mem = 8
    String? modules = "hmmcopy-utils/0.1.1"
  }

  command <<<
    set -o pipefail

    # index
    readCounter \
    --window ~{windowSize} \
    --quality ~{minimumMappingQuality} \
    --chromosome ~{chromosomesToAnalyze} \
    --build ~{bam}

    # convert
    readCounter \
    --window ~{windowSize} \
    --quality ~{minimumMappingQuality} \
    --chromosome ~{chromosomesToAnalyze} \
    ~{bam} | sed "s/chrom=chr/chrom=/" > ~{outputFileNamePrefix}.wig
  >>>

  runtime {
    memory: "~{mem} GB"
    modules: "~{modules}"
  }

  output {
    File wig = "~{outputFileNamePrefix}.wig"
  }
}

# second step : run ichorCNA
task runIchorCNA {
  input{
    File wig
    File? normalWig
    String gcWig = "$ICHORCNA_ROOT/lib/R/ichorCNA/extdata/gc_hg19_1000kb.wig"
    String? mapWig = "$ICHORCNA_ROOT/lib/R/ichorCNA/extdata/map_hg19_1000kb.wig"
    String? normalPanel = "$ICHORCNA_ROOT/lib/R/ichorCNA/extdata/HD_ULP_PoN_1Mb_median_normAutosome_mapScoreFiltered_median.rds"
    String? exonsBed
    String? centromere = "$ICHORCNA_ROOT/lib/R/ichorCNA/extdata/GRCh37.p13_centromere_UCSC-gapTable.txt"
    Float? minMapScore
    Int? rmCentromereFlankLength
    String normal = "\"c(0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9)\""
    String? scStates = "\"c(1, 3)\""
    String? coverage
    String? lambda
    Int? lambdaScaleHyperParam
    String? ploidy = "\"c(2,3)\""
    Int? maxCN = 5
    Boolean? estimateNormal = true
    Boolean? estimateScPrevalence = true
    Boolean? estimatePloidy = true
    Float? maxFracCNASubclone
    Float? maxFracGenomeSubclone
    String? minSegmentBins
    Float? altFracThreshold
    String? chrNormalize
    String? chrTrain = "\"c(1:22)\""
    String? chrs = "\"c(1:22, 'X')\""
    String? genomeBuild
    String? genomeStyle
    Boolean? normalizeMaleX
    Float? fracReadsInChrYForMale
    Boolean? includeHOMD = true
    Float? txnE = 0.9999
    Int? txnStrength = 10000
    String? plotFileType
    String? plotYLim
    String? outDir = "./"
    String? libdir

    String? outputFileNamePrefix = "output"
    String? modules = "ichorcna/0.2"
    Int? mem = 8
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
    ~{"--chrs " + chrs} \
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
}
