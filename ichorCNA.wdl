version 1.0
workflow ichorCNA {

  input {

    String? modules = "ichorcna/0.2 hmmcopy-utils"
    File input_bam
    File input_bam_index
    String sampleName
    Boolean includeHOMD
    Int maxCNstates
    Int window
    Int qual
    Int mem
    String chromosome_wig
  }

call runReadCounter{
  input:
    input_bam=input_bam,
    input_bam_index=input_bam_index,
    sampleName=sampleName,
    window=window,
    qual=qual,
    modules=modules,
    mem=mem,
    chromosome_wig=chromosome_wig
 }

call runIchorCNA {
input:
  includeHOMD=includeHOMD,
  maxCNstates=maxCNstates,
  sampleName=sampleName,
  modules=modules,
  mem=mem,
  bamWig=runReadCounter.bamWig

 }
}


task copyOverBam {
  input {
    File input_bam
    File input_bam_index
  }
  command <<<
  ls -l $(dirname ~{input_bam})
  ls -l $(dirname ~{input_bam_index})
  >>>
}

# First step : read counter
task runReadCounter {
  input{
    File input_bam
    File input_bam_index
    String sampleName
    Int window
    Int qual
    Int mem
    String? modules
    String chromosome_wig
  }
  runtime {
        memory: "~{mem} GB"
        modules: "~{modules}"
    }

  command {
    set -o pipefail
    readCounter \
      --window ~{window} \
      --quality ~{qual} \
      --chromosome "~{chromosome_wig}" \
      -b ~{input_bam}
    $HMMCOPY_UTILS_ROOT/bin/readCounter --window ~{window} \
        --quality ~{qual} \
        --chromosome ~{chromosome_wig} \
        ~{input_bam} | sed "s/chrom=chr/chrom=/" > ~{sampleName}.wig
    }
  output {
    File bamWig = "${sampleName}.wig"
  }
}

# second step : run ichorCNA
task runIchorCNA {
  input{
    File bamWig
    String sampleName
    Boolean includeHOMD
    Int maxCNstates
    String? modules
    Int mem
  }

  runtime {
        memory: "~{mem} GB"
        modules: "~{modules}"
    }

command {

  runIchorCNA \
    --id ~{sampleName} \
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
    --outDir ~{sampleName}_ichorCNA
    }
output {
  File ichorResults = "~{sampleName}_ichorCNA"
  }
}
