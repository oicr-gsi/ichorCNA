version 1.0

workflow ichorCNA {

  input {
    File bam
    File bamIndex
    String? outputFileNamePrefix
    Boolean includeHOMD
    Int maxCNstates
    Int window
    Int qual
    String chromosomeWig
  }

  call runReadCounter{
    input:
      bam=bam,
      bamIndex=bamIndex,
      outputFileNamePrefix=outputFileNamePrefix,
      window=window,
      qual=qual,
      chromosomeWig=chromosomeWig
 }

  call runIchorCNA {
    input:
      includeHOMD=includeHOMD,
      maxCNstates=maxCNstates,
      outputFileNamePrefix=outputFileNamePrefix,
      bamWig=runReadCounter.bamWig
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
    Int window
    Int qual
    Int? mem = 8
    String? modules = "hmmcopy-utils/0.1.1"
    String chromosomeWig
  }

  command <<<
    set -o pipefail

    # index
    readCounter \
    --window ~{window} \
    --quality ~{qual} \
    --chromosome ~{chromosomeWig} \
    --build ~{bam}

    # convert
    readCounter \
    --window ~{window} \
    --quality ~{qual} \
    --chromosome ~{chromosomeWig} \
    ~{bam} | sed "s/chrom=chr/chrom=/" > ~{outputFileNamePrefix}.wig
  >>>

  runtime {
    memory: "~{mem} GB"
    modules: "~{modules}"
  }

  output {
    File bamWig = "~{outputFileNamePrefix}.wig"
  }
}

# second step : run ichorCNA
task runIchorCNA {
  input{
    File bamWig
    String? outputFileNamePrefix = "output"
    Boolean includeHOMD
    Int maxCNstates
    String? modules = "ichorcna/0.2"
    Int? mem = 8
  }

  command <<<
    mkdir ~{outputFileNamePrefix}_ichorCNA

    runIchorCNA \
    --id ~{outputFileNamePrefix} \
    --WIG ~{bamWig} \
    --ploidy "c(2,3)" \
    --normal "c(0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9)" \
    --maxCN 5 \
    --gcWig $ICHORCNA_ROOT/lib/R/ichorCNA/extdata/gc_hg19_1000kb.wig \
    --mapWig $ICHORCNA_ROOT/lib/R/ichorCNA/extdata/map_hg19_1000kb.wig \
    --centromere $ICHORCNA_ROOT/lib/R/ichorCNA/extdata/GRCh37.p13_centromere_UCSC-gapTable.txt \
    --normalPanel $ICHORCNA_ROOT/lib/R/ichorCNA/extdata/HD_ULP_PoN_1Mb_median_normAutosome_mapScoreFiltered_median.rds \
    --includeHOMD ~{true="True" false="False" includeHOMD} \
    --chrs "c(1:22, 'X')" \
    --chrTrain "c(1:22)" \
    --estimateNormal True \
    --estimatePloidy True \
    --estimateScPrevalence True \
    --scStates "c(1, 3)" \
    --txnE 0.9999 \
    --txnStrength 10000 \
    --outDir ./

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
