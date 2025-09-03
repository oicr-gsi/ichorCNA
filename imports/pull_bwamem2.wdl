version 1.0

workflow bwamem2 {
    input {
        File fastqR1
        File? fastqR2
        String outputFileNamePrefix
        Int numChunk = 1
        Boolean doUMIextract = false
        Boolean doTrim = false
        String reference
        Int? numReads
    }

    parameter_meta {
        fastqR1: "Fastq file for read 1"
        fastqR2: "Fastq file for read 2"
        outputFileNamePrefix: "Prefix for output files"
        numChunk: "Number of chunks to split fastq file [1, no splitting]"
        doUMIextract: "If true, UMI will be extracted before alignment [false]"
        doTrim: "If true, adapters will be trimmed before alignment [false]"
        reference: "The genome reference build. For example: hg19, hg38, mm10"
        numReads: "Number of reads"
    }

    Map[String,String] bwamem2_modules_by_genome = { 
    "hg19": "samtools/1.9 bwa-mem2/2.2.1 hg19-bwamem2-index/2.2.1",
    "hg38": "samtools/1.9 bwa-mem2/2.2.1 hg38-bwamem2-index-with-alt/2.2.1"
    }

    Map[String,String] bwamem2Ref_by_genome = { 
    "hg19": "$HG19_BWAMEM2_INDEX_ROOT/hg19_index",
    "hg38": "$HG38_BWAMEM2_INDEX_WITH_ALT_ROOT/hg38_random.fa"
    }


    String bwamem2_modules = bwamem2_modules_by_genome [ reference ]
    String bwamem2_ref = bwamem2Ref_by_genome [ reference ]

    if (numChunk > 1) {
        call countChunkSize {
            input:
            fastqR1 = fastqR1,
            numChunk = numChunk,
            numReads = numReads
        }
    
        call slicer as slicerR1 { 
            input: 
            fastqR = fastqR1,
            chunkSize = countChunkSize.chunkSize
        }
        if (defined(fastqR2)) {
            # workaround for converting File? to File
            File fastqR2_ = select_all([fastqR2])[0]
            call slicer as slicerR2 {
                input:
                fastqR = fastqR2_,
                chunkSize = countChunkSize.chunkSize
            }
        }
    }

    Array[File] fastq1 = select_first([slicerR1.chunkFastq, [fastqR1]])

    if(defined(fastqR2)) {
      Array[File?] fastq2 = select_first([slicerR2.chunkFastq, [fastqR2]])
      Array[Pair[File,File?]] pairedFastqs = zip(fastq1,fastq2)
    }

    if(!defined(fastqR2)) {
      Array[Pair[File,File?]] singleFastqs = cross(fastq1,[fastqR2])
    }

    Array[Pair[File,File?]] outputs = select_first([pairedFastqs, singleFastqs])

    scatter (p in outputs) {

        if (doUMIextract) {
            call extractUMIs { 
                input:
                fastq1 = p.left,
                fastq2 = p.right,
            }

        }

        if (doTrim) {
            call adapterTrimming { 
                input:
                fastqR1 = select_first([extractUMIs.fastqR1, p.left]),
                fastqR2 = if (defined(fastqR2)) then select_first([extractUMIs.fastqR2, p.right]) else fastqR2,
            }
        }


        call runBwamem2 { 
                input: 
                read1s = select_first([adapterTrimming.resultR1, extractUMIs.fastqR1, p.left]),
                read2s = if (defined(fastqR2)) then select_first([adapterTrimming.resultR2, extractUMIs.fastqR2, p.right]) else fastqR2,
                modules = bwamem2_modules,
                bwa2Ref = bwamem2_ref
        }    
    }

    call bamMerge {
        input:
        bams = runBwamem2.outputBam,
        outputFileNamePrefix = outputFileNamePrefix
    }

    call indexBam { 
        input: 
        inputBam = bamMerge.outputMergedBam
    }

    if (doTrim) {
        call adapterTrimmingLog {
            input:
            inputLogs = select_all(adapterTrimming.log),
            outputFileNamePrefix = outputFileNamePrefix,
            numChunk = numChunk,
            singleEnded = if (defined(fastqR2)) then false else true
        }
    }

    meta {
        author: "Xuemei Luo"
        email: "xuemei.luo@oicr.on.ca"
        description: "This workflow aligns sequence data provided as fastq files against a genomic reference using bwamem2.  Prior to alignment, there are options to remove 5' umi sequence and to trim off 3' sequencing adapter. Readgroup information to be injected into the bam header needs to be provided.  The workflow can also split the input data into a requested number of chunks, align each separately then merge the separate alignments into a single bam file.  This decreases the workflow run time.  Optional bwa mem parameters can be provided to the workflow."
        dependencies: [
        {
            name: "bwa-mem2/2.2.1",
            url: "https://github.com/bwa-mem2/bwa-mem2/releases/download/v2.2.1/bwa-mem2-2.2.1_x64-linux.tar.bz2"
        },
        {
            name: "samtools/1.9",
            url: "https://github.com/samtools/samtools/archive/0.1.19.tar.gz"
        },
        {
            name: "cutadapt/1.8.3",
            url: "https://cutadapt.readthedocs.io/en/v1.8.3/"
        },
        {
            name: "slicer/0.3.0",
            url: "https://github.com/OpenGene/slicer/archive/v0.3.0.tar.gz"
        },
        {
            name: "python/3.7",
            url: "https://www.python.org"
        },      
        {
            name: "barcodex-rs/0.1.2",
            url: "https://github.com/oicr-gsi/barcodex-rs/archive/v0.1.2.tar.gz"
        },
        {
            name: "rust/1.2",
            url: "https://www.rust-lang.org/tools/install"
        },
        { 
          name: "gsi software modules : samtools/1.9 bwa-mem2/2.2.1",
          url: "https://gitlab.oicr.on.ca/ResearchIT/modulator"
        },
        { 
          name: "gsi hg38 modules : hg38-bwamem2-index-with-alt/2.2.1",
          url: "https://gitlab.oicr.on.ca/ResearchIT/modulator"
        },
        {
          name: "gsi hg19 modules : hg19-bwamem2-index/2.2.1",
          url: "https://gitlab.oicr.on.ca/ResearchIT/modulator"
        }
      ]
      output_meta: {
       bwamem2Bam: {
          description: "Output Alignment BAM file",
          vidarr_label: "bwamem2Bam"
       },
       bwamem2Index: {
          description: "Index of the Output Alignment file, BAI",
          vidarr_label: "bwamem2Index"
       },
       log: {
          description: "Optional log file",
          vidarr_label: "log"
       },
       cutAdaptAllLogs: {
          description: "Optional log file for cutAdapt",
          vidarr_label: "cutAdaptAllLogs"
       } 
      }
    }

    output {
        File bwamem2Bam = bamMerge.outputMergedBam
        File bwamem2Index = indexBam.outputBai
        File? log = adapterTrimmingLog.summaryLog
        File? cutAdaptAllLogs = adapterTrimmingLog.allLogs
    }
}


task countChunkSize{
    input {
        File fastqR1
        Int numChunk
        Int? numReads
        String modules = "python/3.7"
        Int jobMemory = 16
        Int timeout = 48
    }
    
    parameter_meta {
        fastqR1: "Fastq file for read 1"
        numChunk: "Number of chunks to split fastq file"
        numReads: "Number of reads"
        modules: "Required environment modules"
        jobMemory: "Memory allocated for this job"
        timeout: "Hours before task timeout"
    }
    
    command <<<
        set -euo pipefail

        if [ -z "~{numReads}" ]; then
            totalLines=$(zcat ~{fastqR1} | wc -l)
        else totalLines=$((~{numReads}*4))
        fi
        
        python3 -c "from math import ceil; print (int(ceil(($totalLines/4.0)/~{numChunk})*4))"
    >>>
    
    runtime {
        memory: "~{jobMemory} GB"
        modules: "~{modules}"
        timeout: "~{timeout}"
    }
    
    output {
        String chunkSize =read_string(stdout())
    }

    meta {
    output_meta: {
    chunkSize: "output number of lines per chunk"
    }
    }    
   
}

task slicer {
    input {
        File fastqR         
        String chunkSize
        String modules = "slicer/0.3.0"
        Int jobMemory = 16
        Int timeout = 48
    }
    
    parameter_meta {
        fastqR: "Fastq file"
        chunkSize: "Number of lines per chunk"
        modules: "Required environment modules"
        jobMemory: "Memory allocated for this job"
        timeout: "Hours before task timeout"
    }
    
    command <<<
        set -euo pipefail
        slicer -i ~{fastqR} -l ~{chunkSize} --gzip 
    >>>
    
    runtime {
        memory: "~{jobMemory} GB"
        modules: "~{modules}"
        timeout: "~{timeout}"
    } 
    
    output {
        Array[File] chunkFastq = glob("*.fastq.gz")
    }

    meta {
        output_meta: {
            chunkFastq: "output fastq chunks"
        }
    } 
  
}


task extractUMIs {
        input {
            String umiList = "umiList"
            String outputPrefix = "extractUMIs_output"
            File fastq1
            File? fastq2
            String pattern1 = "(?P<umi_1>^[ACGT]{3}[ACG])(?P<discard_1>T)|(?P<umi_2>^[ACGT]{3})(?P<discard_2>T)"
            String pattern2 = "(?P<umi_1>^[ACGT]{3}[ACG])(?P<discard_1>T)|(?P<umi_2>^[ACGT]{3})(?P<discard_2>T)"
            String modules = "barcodex-rs/0.1.2 rust/1.45.1"
            Int jobMemory = 24
            Int timeout = 12
        }

        parameter_meta {
            umiList: "Reference file with valid UMIs"
            outputPrefix: "Specifies the start of the output files"
            fastq1: "FASTQ file containing read 1"
            fastq2: "FASTQ file containing read 2"
            pattern1: "UMI RegEx pattern 1"
            pattern2: "UMI RegEx pattern 2"
            modules: "Required environment modules"
            jobMemory: "Memory allocated for this job"
            timeout: "Time in hours before task timeout"
        }

        command <<<
            set -euo pipefail

            barcodex-rs --umilist ~{umiList} --prefix ~{outputPrefix} --separator "__" inline \
            --pattern1 '~{pattern1}' --r1-in ~{fastq1} \
            ~{if (defined(fastq2)) then "--pattern2 '~{pattern2}' --r2-in ~{fastq2} " else ""}

            cat ~{outputPrefix}_UMI_counts.json > umiCounts.txt

            tr [,] ',\n' < umiCounts.txt | sed 's/[{}]//' > tmp.txt
            echo "{$(sort -i tmp.txt)}" > new.txt
            tr '\n' ',' < new.txt | sed 's/,$//' > ~{outputPrefix}_UMI_counts.json
        >>>

        runtime {
            modules: "~{modules}"
            memory: "~{jobMemory} GB"
            timeout: "~{timeout}"
        }

        output {
            File fastqR1 = "~{outputPrefix}_R1.fastq.gz"
            File? fastqR2 = "~{outputPrefix}_R2.fastq.gz"
            File discardR1 = "~{outputPrefix}_R1.discarded.fastq.gz"
            File? discardR2 = "~{outputPrefix}_R2.discarded.fastq.gz"
            File extractR1 = "~{outputPrefix}_R1.extracted.fastq.gz"
            File? extractR2 = "~{outputPrefix}_R2.extracted.fastq.gz"
            File umiCounts = "~{outputPrefix}_UMI_counts.json"
            File extractionMetrics = "~{outputPrefix}_extraction_metrics.json"
        }

        meta {
            output_meta: {
                fastqR1: "Read 1 fastq file with UMIs extracted",
                fastqR2: "Read 2 fastq file with UMIs extracted",
                discardR1: "Reads without a matching UMI pattern in read 1",
                discardR2: "Reads without a matching UMI pattern in read 2",
                extractR1: "Extracted reads (UMIs and any spacer sequences) from read 1",
                extractR2: "Extracted reads (UMIs and any spacer sequences) from read 2",
                umiCounts: "Record of UMI counts after extraction",
                extractionMetrics: "Metrics relating to extraction process"
            }
        }
}


task adapterTrimming {
    input {
        File fastqR1
        File? fastqR2
        String modules = "cutadapt/1.8.3"
        Boolean doUMItrim = false
        Int umiLength = 5
        Int trimMinLength = 1
        Int trimMinQuality = 0
        String adapter1 = "AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC"
        String adapter2 = "AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT" 
        String? addParam
        Int jobMemory = 16
        Int timeout = 48  
    }
    
    parameter_meta {
        fastqR1: "Fastq file for read 1"
        fastqR2: "Fastq file for read 2"
        doUMItrim: "If true, do umi trimming"
        umiLength: "The number of bases to trim when doUMItrim is true. If the given length is positive, the bases are removed from the beginning of each read. If it is negative, the bases are removed from the end"
        trimMinLength: "Minimum length of reads to keep"
        trimMinQuality: "Minimum quality of read ends to keep"
        adapter1: "Adapter sequence to trim from read 1"
        adapter2: "Adapter sequence to trim from read 2"
        modules: "Required environment modules"
        addParam: "Additional cutadapt parameters"
        jobMemory: "Memory allocated for this job"
        timeout: "Hours before task timeout"
    }
   
    Array[File] inputs = select_all([fastqR1,fastqR2])
    String resultFastqR1 = "~{basename(fastqR1, ".fastq.gz")}.trim.fastq.gz"
    String resultFastqR2 = if (length(inputs) > 1) then "~{basename(inputs[1], ".fastq.gz")}.trim.fastq.gz" else "None"
    String resultLog = "~{basename(fastqR1, ".fastq.gz")}.log"
    
    command <<<
        set -euo pipefail

        cutadapt -q ~{trimMinQuality} \
                -m ~{trimMinLength} \
                -a ~{adapter1} \
                -o ~{resultFastqR1} \
                ~{if (defined(fastqR2)) then "-A ~{adapter2} -p ~{resultFastqR2} " else ""} \
                ~{if (doUMItrim) then "-u ~{umiLength} -U ~{umiLength} " else ""} \
                ~{addParam} \
                ~{fastqR1} \
                ~{fastqR2} > ~{resultLog}

    >>>
    
    runtime {
        memory: "~{jobMemory} GB"
        modules: "~{modules}"
        timeout: "~{timeout}"
    } 
    
    output { 
        File resultR1 = "~{resultFastqR1}"
        File? resultR2 = "~{resultFastqR2}"
        File log =  "~{resultLog}"     
    }

    meta {
        output_meta: {
            resultR1: "output fastq read 1 after trimming",
            resultR2: "output fastq read 2 after trimming",
            log: "output adpater trimming log"
        }
    } 
}


task runBwamem2 {
    input {
        File read1s
        File? read2s
        String readGroups
        String modules
        String bwa2Ref
        String? addParam
        Int threads = 8
        Int jobMemory = 32
        Int timeout = 96
    }

    parameter_meta {
        read1s: "Fastq file for read 1"
        read2s: "Fastq file for read 2"
        readGroups: "The readgroup information to be injected into the bam header"
        bwa2Ref: "The reference genome to align the sample with by bwa-mem2"
        modules: "Required environment modules"
        addParam: "Additional BWA parameters"
        threads: "Requested CPU threads"
        jobMemory: "Memory allocated for this job"
        timeout: "Hours before task timeout"
    }
    
    String resultBam = "~{basename(read1s)}.bam"
    String tmpDir = "tmp/"

    command <<<
        set -euo pipefail
        mkdir -p ~{tmpDir}
        bwa-mem2 mem -M \
            -t ~{threads} ~{addParam}  \
            -R  ~{readGroups} \
            ~{bwa2Ref} \
            ~{read1s} \
            ~{read2s} \
        | \
        samtools sort -O bam -T ~{tmpDir} -o ~{resultBam} - 
    >>>

    runtime {
        modules: "~{modules}"
        memory:  "~{jobMemory} GB"
        cpu:     "~{threads}"
        timeout: "~{timeout}"
    }  
    
    output {
        File outputBam = "~{resultBam}"
    }

    meta {
        output_meta: {
            outputBam: "output bam aligned to genome"
        }
    }

}

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

task indexBam {
    input {
        File  inputBam
        Int   jobMemory = 12
        String modules  = "samtools/1.9"
        Int timeout     = 48
    }
    parameter_meta {
        inputBam:  "Input bam file"
        jobMemory: "Memory allocated indexing job"
        modules:   "Modules for running indexing job"
        timeout:   "Hours before task timeout"
    }

    String resultBai = "~{basename(inputBam)}.bai"

    command <<<
        set -euo pipefail
        samtools index ~{inputBam} ~{resultBai}
    >>>

    runtime {
        memory: "~{jobMemory} GB"
        modules: "~{modules}"
        timeout: "~{timeout}"
    }

    output {
        File outputBai = "~{resultBai}"
    }

    meta {
        output_meta: {
            outputBai: "output index file for bam aligned to genome"
        }
    }

}

task adapterTrimmingLog {
    input {
        Array[File] inputLogs
        String outputFileNamePrefix
        Int   numChunk
        Boolean singleEnded = false
        Int   jobMemory = 12
        Int timeout     = 48


    }
    parameter_meta {
        inputLogs:  "Input log files"
        outputFileNamePrefix: "Prefix for output file"
        numChunk: "Number of chunks to split fastq file"
        singleEnded: "true if reads are single ended"
        jobMemory: "Memory allocated indexing job"
        timeout:   "Hours before task timeout"
    }

    String allLog = "~{outputFileNamePrefix}.txt"
    String log = "~{outputFileNamePrefix}.log"

    command <<<
        set -euo pipefail
        awk 'BEGINFILE {print "###################################\n"}{print}' ~{sep=" " inputLogs} > ~{allLog}

        totalBP=$(cat ~{allLog} | grep "Total basepairs processed:" | cut -d":" -f2 | sed 's/^[ \t]*//; s/ bp//; s/,//g' | awk '{x+=$1}END{print x}')

        bpQualitytrimmed=$(cat ~{allLog} | grep "Quality-trimmed:" | cut -d":" -f2 | sed 's/^[ \t]*//; s/ bp//; s/,//g; s/ (.*)//' | awk '{x+=$1}END{print x}')
        percentQualitytrimmed=$(awk -v A="${bpQualitytrimmed}" -v B="${totalBP}" 'BEGIN {printf "%0.1f\n", A*100.0/B}')

        bpTotalWritten=$(cat ~{allLog} | grep "Total written (filtered):" | cut -d":" -f2 | sed 's/^[ \t]*//; s/ bp//; s/,//g; s/ (.*)//' | awk '{x+=$1}END{print x}')
        percentBPWritten=$(awk -v A="${bpTotalWritten}" -v B="${totalBP}" 'BEGIN {printf "%0.1f\n", A*100.0/B}')

        echo -e "This is a cutadapt summary from ~{numChunk} fastq chunks\n" > ~{log}

        if ! ~{singleEnded} ; then
          totalRead=$(cat ~{allLog} | grep "Total read pairs processed:" | cut -d":" -f2 | sed 's/ //g; s/,//g' | awk '{x+=$1}END{print x}')
          adapterR1=$(cat ~{allLog} | grep " Read 1 with adapter:" | cut -d ":" -f2 | sed 's/^[ \t]*//; s/ (.*)//; s/,//g'| awk '{x+=$1}END{print x}')
          percentAdapterR1=$(awk -v A="${adapterR1}" -v B="${totalRead}" 'BEGIN {printf "%0.1f\n", A*100.0/B}')
          adapterR2=$(cat ~{allLog} | grep " Read 2 with adapter:" | cut -d ":" -f2 | sed 's/^[ \t]*//; s/ (.*)//; s/,//g'| awk '{x+=$1}END{print x}')
          percentAdapterR2=$(awk -v A="${adapterR2}" -v B="${totalRead}" 'BEGIN {printf "%0.1f\n", A*100.0/B}')

          shortPairs=$(cat ~{allLog} | grep "Pairs that were too short:" | cut -d ":" -f2 | sed 's/^[ \t]*//; s/ (.*)//; s/,//g'| awk '{x+=$1}END{print x}')
          percentShortPairs=$(awk -v A="${shortPairs}" -v B="${totalRead}" 'BEGIN {printf "%0.1f\n", A*100.0/B}')

          pairsWritten=$(cat ~{allLog} | grep "Pairs written (passing filters): " | cut -d ":" -f2 | sed 's/^[ \t]*//; s/ (.*)//; s/,//g'| awk '{x+=$1}END{print x}')
          percentpairsWritten=$(awk -v A="${pairsWritten}" -v B="${totalRead}" 'BEGIN {printf "%0.1f\n", A*100.0/B}')

          bpR1=$(cat ~{allLog} | grep -A 2 "Total basepairs processed:" | grep "Read 1:" | cut -d":" -f2 | sed 's/^[ \t]*//; s/ bp//; s/,//g' | awk '{x+=$1}END{print x}')
          bpR2=$(cat ~{allLog} | grep -A 2 "Total basepairs processed:" | grep "Read 2:" | cut -d":" -f2 | sed 's/^[ \t]*//; s/ bp//; s/,//g' | awk '{x+=$1}END{print x}')

          bpQualitytrimmedR1=$(cat ~{allLog} | grep -A 2 "Quality-trimmed:" | grep "Read 1:" | cut -d":" -f2 | sed 's/^[ \t]*//; s/ bp//; s/,//g' | awk '{x+=$1}END{print x}')
          bpQualitytrimmedR2=$(cat ~{allLog} | grep -A 2 "Quality-trimmed:" | grep "Read 2:" | cut -d":" -f2 | sed 's/^[ \t]*//; s/ bp//; s/,//g' | awk '{x+=$1}END{print x}')

          bpWrittenR1=$(cat ~{allLog} | grep -A 2 "Total written (filtered):" | grep "Read 1:" | cut -d":" -f2 | sed 's/^[ \t]*//; s/ bp//; s/,//g' | awk '{x+=$1}END{print x}')
          bpWrittenR2=$(cat ~{allLog} | grep -A 2 "Total written (filtered):" | grep "Read 2:" | cut -d":" -f2 | sed 's/^[ \t]*//; s/ bp//; s/,//g' | awk '{x+=$1}END{print x}')

          echo -e "Total read pairs processed:\t${totalRead}" >> ~{log}
          echo -e "  Read 1 with adapter:\t${adapterR1} (${percentAdapterR1}%)" >> ~{log}
          echo -e "  Read 2 with adapter:\t${adapterR2} (${percentAdapterR2}%)" >> ~{log}
          echo -e "Pairs that were too short:\t${shortPairs} (${percentShortPairs}%)" >> ~{log}
          echo -e "Pairs written (passing filters):\t${pairsWritten} (${percentpairsWritten}%)\n\n" >> ~{log}
          echo -e "Total basepairs processed:\t${totalBP} bp" >> ~{log}
          echo -e "  Read 1:\t${bpR1} bp" >> ~{log}
          echo -e "  Read 2:\t${bpR2} bp" >> ~{log}
          echo -e "Quality-trimmed:\t${bpQualitytrimmed} bp (${percentQualitytrimmed}%)" >> ~{log}
          echo -e "  Read 1:\t${bpQualitytrimmedR1} bp" >> ~{log}
          echo -e "  Read 2:\t${bpQualitytrimmedR2} bp" >> ~{log}
          echo -e "Total written (filtered):\t${bpTotalWritten} bp (${percentBPWritten}%)" >> ~{log}
          echo -e "  Read 1:\t${bpWrittenR1} bp" >> ~{log}
          echo -e "  Read 2:\t${bpWrittenR2} bp" >> ~{log}

        else 
          totalRead=$(cat ~{allLog} | grep "Total reads processed:" | cut -d":" -f2 | sed 's/ //g; s/,//g' | awk '{x+=$1}END{print x}')
          adapterR=$(cat ~{allLog} | grep "Reads with adapters:" | cut -d ":" -f2 | sed 's/^[ \t]*//; s/ (.*)//; s/,//g'| awk '{x+=$1}END{print x}')
          percentAdapterR=$(awk -v A="${adapterR}" -v B="${totalRead}" 'BEGIN {printf "%0.1f\n", A*100.0/B}')

          shortReads=$(cat ~{allLog} | grep "Reads that were too short:" | cut -d ":" -f2 | sed 's/^[ \t]*//; s/ (.*)//; s/,//g'| awk '{x+=$1}END{print x}')
          percentShortReads=$(awk -v A="${shortReads}" -v B="${totalRead}" 'BEGIN {printf "%0.1f\n", A*100.0/B}')

          ReadsWritten=$(cat ~{allLog} | grep "Reads written (passing filters): " | cut -d ":" -f2 | sed 's/^[ \t]*//; s/ (.*)//; s/,//g'| awk '{x+=$1}END{print x}')
          percentreadsWritten=$(awk -v A="${ReadsWritten}" -v B="${totalRead}" 'BEGIN {printf "%0.1f\n", A*100.0/B}')                 

          echo -e "Total reads processed:\t${totalRead}" >> ~{log}
          echo -e "Reads with adapters:\t${adapterR} (${percentAdapterR}%)" >> ~{log}
          echo -e "Reads that were too short:\t${shortReads} (${percentShortReads}%)" >> ~{log}
          echo -e "Reads written (passing filters):\t${ReadsWritten} (${percentreadsWritten}%)\n\n" >> ~{log}
          echo -e "Total basepairs processed:\t${totalBP} bp" >> ~{log}
          echo -e "Quality-trimmed:\t${bpQualitytrimmed} bp (${percentQualitytrimmed}%)" >> ~{log}
          echo -e "Total written (filtered):\t${bpTotalWritten} bp (${percentBPWritten}%)" >> ~{log}
        fi
    >>>

    runtime {
        memory: "~{jobMemory} GB"
        timeout: "~{timeout}"
    }
  
    output {
        File summaryLog = "~{log}"
        File allLogs = "~{allLog}"
    }

    meta {
        output_meta: {
            summaryLog: "a summary log file for adapter trimming",
            allLogs: "a file containing all logs for adapter trimming for each fastq chunk"
        }
    }

}
