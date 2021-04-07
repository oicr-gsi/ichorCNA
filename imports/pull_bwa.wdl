version 1.0

workflow bwaMem {
    input {
        Int adapterTrimmingLog_timeout = 48
        Int adapterTrimmingLog_jobMemory = 12
        Int indexBam_timeout = 48
        String indexBam_modules = "samtools/1.9"
        Int indexBam_jobMemory = 12
        Int bamMerge_timeout = 72
        String bamMerge_modules = "samtools/1.9"
        Int bamMerge_jobMemory = 32
        Int runBwaMem_timeout = 96
        Int runBwaMem_jobMemory = 32
        Int runBwaMem_threads = 8
        String? runBwaMem_addParam
        String runBwaMem_bwaRef
        String runBwaMem_modules
        Int adapterTrimming_timeout = 48
        Int adapterTrimming_jobMemory = 16
        String? adapterTrimming_addParam
        String adapterTrimming_modules = "cutadapt/1.8.3"
        Int slicerR2_timeout = 48
        Int slicerR2_jobMemory = 16
        String slicerR2_modules = "slicer/0.3.0"
        Int slicerR1_timeout = 48
        Int slicerR1_jobMemory = 16
        String slicerR1_modules = "slicer/0.3.0"
        Int countChunkSize_timeout = 48
        Int countChunkSize_jobMemory = 16
        File fastqR1
        File? fastqR2
        String readGroups
        String outputFileNamePrefix = "output"
        Int numChunk = 1
        Boolean doTrim = false
        Int trimMinLength = 1
        Int trimMinQuality = 0
        String adapter1 = "AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC"
        String adapter2 = "AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT"
    }

    parameter_meta {
        adapterTrimmingLog_timeout: "Hours before task timeout"
        adapterTrimmingLog_jobMemory: "Memory allocated indexing job"
        indexBam_timeout: "Hours before task timeout"
        indexBam_modules: "Modules for running indexing job"
        indexBam_jobMemory: "Memory allocated indexing job"
        bamMerge_timeout: "Hours before task timeout"
        bamMerge_modules: "Required environment modules"
        bamMerge_jobMemory: "Memory allocated indexing job"
        runBwaMem_timeout: "Hours before task timeout"
        runBwaMem_jobMemory: "Memory allocated for this job"
        runBwaMem_threads: "Requested CPU threads"
        runBwaMem_addParam: "Additional BWA parameters"
        runBwaMem_bwaRef: "The reference genome to align the sample with by BWA"
        runBwaMem_modules: "Required environment modules"
        adapterTrimming_timeout: "Hours before task timeout"
        adapterTrimming_jobMemory: "Memory allocated for this job"
        adapterTrimming_addParam: "Additional cutadapt parameters"
        adapterTrimming_modules: "Required environment modules"
        slicerR2_timeout: "Hours before task timeout"
        slicerR2_jobMemory: "Memory allocated for this job"
        slicerR2_modules: "Required environment modules"
        slicerR1_timeout: "Hours before task timeout"
        slicerR1_jobMemory: "Memory allocated for this job"
        slicerR1_modules: "Required environment modules"
        countChunkSize_timeout: "Hours before task timeout"
        countChunkSize_jobMemory: "Memory allocated for this job"
        fastqR1: "fastq file for read 1"
        fastqR2: "fastq file for read 2"
        readGroups: "Complete read group header line"
        outputFileNamePrefix: "Prefix for output file"
        numChunk: "number of chunks to split fastq file [1, no splitting]"
        doTrim: "if true, adapters will be trimmed before alignment"
        trimMinLength: "minimum length of reads to keep [1]"
        trimMinQuality: "minimum quality of read ends to keep [0]"
        adapter1: "adapter sequence to trim from read 1 [AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC]"
        adapter2: "adapter sequence to trim from read 2 [AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT]"

    }

    if (numChunk > 1) {
        call countChunkSize {
            input:
            timeout = countChunkSize_timeout,
            jobMemory = countChunkSize_jobMemory,
            fastqR1 = fastqR1,
            numChunk = numChunk
        }
    
        call slicer as slicerR1 { 
            input: 
            timeout = slicerR1_timeout,
            jobMemory = slicerR1_jobMemory,
            modules = slicerR1_modules,
            fastqR = fastqR1,
            chunkSize = countChunkSize.chunkSize
        }
        if (defined(fastqR2)) {
            # workaround for converting File? to File
            File fastqR2_ = select_all([fastqR2])[0]
            call slicer as slicerR2 {
                input:
                timeout = slicerR2_timeout,
                jobMemory = slicerR2_jobMemory,
                modules = slicerR2_modules,
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
        if (doTrim) {
            call adapterTrimming { 
                input:
                timeout = adapterTrimming_timeout,
                jobMemory = adapterTrimming_jobMemory,
                addParam = adapterTrimming_addParam,
                modules = adapterTrimming_modules,
                fastqR1 = p.left,
                fastqR2 = p.right,
                trimMinLength = trimMinLength,
                trimMinQuality = trimMinQuality,
                adapter1 = adapter1,
                adapter2 = adapter2
            }
        }
        call runBwaMem  { 
                input: 
                timeout = runBwaMem_timeout,
                jobMemory = runBwaMem_jobMemory,
                threads = runBwaMem_threads,
                addParam = runBwaMem_addParam,
                bwaRef = runBwaMem_bwaRef,
                modules = runBwaMem_modules,
                read1s = select_first([adapterTrimming.resultR1, p.left]),
                read2s = if (defined(fastqR2)) then select_first([adapterTrimming.resultR2, p.right]) else fastqR2,
                readGroups = readGroups
        }    
    }

    call bamMerge {
        input:
        timeout = bamMerge_timeout,
        modules = bamMerge_modules,
        jobMemory = bamMerge_jobMemory,
        bams = runBwaMem.outputBam,
        outputFileNamePrefix = outputFileNamePrefix
    }

    call indexBam { 
        input: 
        timeout = indexBam_timeout,
        modules = indexBam_modules,
        jobMemory = indexBam_jobMemory,
        inputBam = bamMerge.outputMergedBam
    }

    if (doTrim) {
        call adapterTrimmingLog {
            input:
            timeout = adapterTrimmingLog_timeout,
            jobMemory = adapterTrimmingLog_jobMemory,
            inputLogs = select_all(adapterTrimming.log),
            outputFileNamePrefix = outputFileNamePrefix,
            numChunk = numChunk,
            singleEnded = if (defined(fastqR2)) then false else true
        }
    }

    meta {
        author: "Xuemei Luo"
        email: "xuemei.luo@oicr.on.ca"
        description: "BwaMem Workflow version 2.0"
        dependencies: [
        {
            name: "bwa/0.7.12",
            url: "https://github.com/lh3/bwa/archive/0.7.12.tar.gz"
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
        }
      ]
    }

    output {
        File bwaMemBam = bamMerge.outputMergedBam
        File bwaMemIndex = indexBam.outputBai
        File? log = adapterTrimmingLog.summaryLog
        File? cutAdaptAllLogs = adapterTrimmingLog.allLogs
    }
}


task countChunkSize{
    input {
        File fastqR1
        Int numChunk
        Int jobMemory = 16
        Int timeout = 48
    }
    
    parameter_meta {
        fastqR1: "Fastq file for read 1"
        numChunk: "Number of chunks to split fastq file"
        jobMemory: "Memory allocated for this job"
        timeout: "Hours before task timeout"
    }
    
    command <<<
        set -euo pipefail
        totalLines=$(zcat ~{fastqR1} | wc -l)
        python -c "from math import ceil; print int(ceil(($totalLines/4.0)/~{numChunk})*4)"
    >>>
    
    runtime {
        memory: "~{jobMemory} GB"
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

task adapterTrimming {
    input {
        File fastqR1
        File? fastqR2
        String modules = "cutadapt/1.8.3"
        Int trimMinLength
        Int trimMinQuality
        String adapter1
        String adapter2
        String? addParam
        Int jobMemory = 16
        Int timeout = 48
    }
    
    parameter_meta {
        fastqR1: "Fastq file for read 1"
        fastqR2: "Fastq file for read 2"
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
            -a ~{adapter1}  \
            -o ~{resultFastqR1} \
            ~{if (defined(fastqR2)) then "-A ~{adapter2} -p ~{resultFastqR2} " else ""} \
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


task runBwaMem {
    input {
        File read1s
        File? read2s
        String readGroups
        String modules
        String bwaRef
        String? addParam
        Int threads = 8
        Int jobMemory = 32
        Int timeout = 96
    }

    parameter_meta {
        read1s: "Fastq file for read 1"
        read2s: "Fastq file for read 2"
        readGroups: "Array of readgroup lines"
        bwaRef: "The reference genome to align the sample with by BWA"
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
        bwa mem -M \
            -t ~{threads} ~{addParam}  \
            -R  ~{readGroups} \
            ~{bwaRef} \
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
 

