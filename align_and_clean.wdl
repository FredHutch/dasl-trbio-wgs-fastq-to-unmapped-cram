version 1.0
## WGS data preprocessing workflow for downstream variant calling.
## Input requirements:
## - Pair-end sequencing data in unmapped BAM (uBAM) format that comply with the following requirements:
## - - filenames all have the same suffix (we use ".unmapped.bam")
## - - files must pass validation by ValidateSamFile (a Picard tool)
## - - reads are provided in query-sorted order
## - - all reads must have an RG tag
##
## Output Files:
## - An analysis-ready recalibrated bam and it's index
## - QC stats from Picard 
## 
workflow WGS_preprocess_for_variants {
  input {

    File? pon
    File? pon_index
    # Reference Genome information
    String ref_name
    File ref_fasta
    File ref_fasta_index
    File ref_dict
    # This is the .alt file from bwa-kit (https://github.com/lh3/bwa/tree/master/bwakit),
    # listing the reference contigs that are "alternative". Leave blank in JSON for legacy
    # references such as b37 and hg19.
    File? ref_alt
    File ref_amb
    File ref_ann
    File ref_bwt
    File ref_pac
    File ref_sa
    File dbSNP_vcf
    File dbSNP_vcf_index
    Array[File] known_indels_sites_VCFs
    Array[File] known_indels_sites_indices
  }

    # Sample information
    # These aren't currently used but it might be a good idea to do so
    #String sampleName = "refForTest"
    #String sampleID = "reftest"
    # for unique IDs for datasets, if used during unmapped bam workflow
    #File bamFile = "s3://fh-pi-paguirigan-a/tg/TGR-Analyses/PairedFastqsToUnmappedBam/b2e28294-440a-412f-ad1b-ca29d20f06b9/1690-98-BULK_M00000865.unmapped.bam"
    File bamFile = "s3://fh-pi-paguirigan-a/tg/TGR-Analyses/PairedFastqsToUnmappedBam/b2e28294-440a-412f-ad1b-ca29d20f06b9/1690-98-MYE_M00000866.unmapped.bam"


    # Docker containers this workflow has been designed for
    String GATKdocker = "broadinstitute/gatk:4.1.8.0"
    String bwadocker = "fredhutch/bwa:0.7.17"

    File wgs_intervals = "s3://fh-ctr-public-reference-data/genome_data/human/hg38/resources_broad_hg38_v0_wgs_calling_regions.hg38.interval_list"
    Int scatter_count = 30

    # Get the basename, i.e. strip the filepath and the extension
    String bam_basename = basename(bamFile, ".unmapped.bam")
    String base_file_name = bam_basename + "." + ref_name

    Int bwaThreads = 16

  # Convert unmapped bam to interleaved fastq
  call splitCram {
    input:
      input_bam = bamFile,
      base_file_name = base_file_name,
      taskDocker = GATKdocker
  }

  #  Map reads to reference
  call cramToBwaMemMerge {
    input:
      input_fastq = SamToFastq.output_fastq,
      base_file_name = base_file_name,
      ref_fasta = ref_fasta,
      ref_fasta_index = ref_fasta_index,
      ref_dict = ref_dict,
      ref_alt = ref_alt,
      ref_amb = ref_amb,
      ref_ann = ref_ann,
      ref_bwt = ref_bwt,
      ref_pac = ref_pac,
      ref_sa = ref_sa,
      bwaThreads = bwaThreads,
      taskDocker = bwadocker
  }

    # Aggregate aligned+merged flowcell BAM files and mark duplicates
  call MarkDuplicatesSpark {
    input:
      input_bam = MergeBamAlignment.output_bam,
      output_bam_basename = base_file_name + ".aligned.duplicates_marked",
      metrics_filename = base_file_name + ".duplicate_metrics",
      taskDocker = GATKdocker
  }

  call SplitIntervals {
    input:
      intervals = wgs_intervals,
      ref_fasta = ref_fasta,
      ref_fasta_index = ref_fasta_index,
      ref_dict = ref_dict,
      scatter_count = scatter_count,
      taskDocker = GATKdocker
  }


  # Perform Base Quality Score Recalibration (BQSR) on the sorted BAM in parallel
  scatter (subinterval in SplitIntervals.interval_files) {

  # Generate the recalibration model by interval and apply it
  call ApplyBaseRecalibrator {
    input:
      input_bam = MarkDuplicatesSpark.output_bam,
      input_bam_index = MarkDuplicatesSpark.output_bai,
      base_file_name = base_file_name,
      dbSNP_vcf = dbSNP_vcf,
      dbSNP_vcf_index = dbSNP_vcf_index,
      known_indels_sites_VCFs = known_indels_sites_VCFs,
      known_indels_sites_indices = known_indels_sites_indices,
      ref_dict = ref_dict,
      ref_fasta = ref_fasta,
      ref_fasta_index = ref_fasta_index,
      intervals = subinterval,
      taskDocker = GATKdocker
    }
  }
  # Merge the recalibrated BAM files resulting from by-interval recalibration
  call GatherBamFiles {
    input:
      input_bams = ApplyBaseRecalibrator.recalibrated_bam,
      base_file_name = base_file_name,
      taskDocker = GATKdocker
  }
    # call CollectAlignmentSummaryMetrics {
    #   input: 
    #     input_bam = sampleApplyBaseRecalibrator.recalibrated_bam,
    #     base_file_name = base_file_name,
    #     ref_fasta = ref_fasta,
    #     ref_fasta_index = ref_fasta_index,
    #     docker = GATKdocker
    # }
    # call CollectWgsMetrics {
    #   input: 
    #     input_bam = sampleApplyBaseRecalibrator.recalibrated_bam,
    #     base_file_name = base_file_name,
    #     ref_fasta = ref_fasta,
    #     ref_fasta_index = ref_fasta_index,
    #     docker = GATKdocker
    # }

  # Outputs that will be retained when execution is complete
  output {
    File analysisReadyBam = GatherBamFiles.output_bam
    File analysisReadyIndex = GatherBamFiles.output_bai
    #File wgsMetrics = CollectWgsMetrics.out
    #File alignmentMetrics = CollectAlignmentSummaryMetrics.out
  }
# End workflow
}
#### TASK DEFINITIONS

# Combine multiple recalibrated BAM files from scattered ApplyRecalibration runs
task GatherBamFiles {
  input {
    Array[File] input_bams
    String base_file_name
    String taskDocker
  }

  command {
    gatk --java-options "-Dsamjdk.compression_level=5 -Xms3g" \
      GatherBamFiles \
      --INPUT ~{sep=' --INPUT ' input_bams} \
      --OUTPUT ~{base_file_name}.bam \
      --CREATE_INDEX true \
      --CREATE_MD5_FILE true
  }
  runtime {
    docker: taskDocker
    memory: "4 GB"
    cpu: 2
  }
  output {
    File output_bam = "~{base_file_name}.bam"
    File output_bai = "~{base_file_name}.bai"
    File output_bam_md5 = "~{base_file_name}.bam.md5"
  }
}
# Read unmapped BAM, convert to FASTQ
task SamToFastq {
  input {
    File input_bam
    String base_file_name
    String taskDocker
  }
  command {
    set -eo pipefail

    gatk --java-options "-Dsamjdk.compression_level=5 -Xms8g" \
      SamToFastq \
      --INPUT ~{input_bam} \
      --FASTQ ~{base_file_name}.fastq.gz \
      --INTERLEAVE true \
      --INCLUDE_NON_PF_READS true 
  }
  output {
    File output_fastq = "~{base_file_name}.fastq.gz"
  }
  runtime {
    memory: "16 GB"
    docker: taskDocker
    walltime: "2:00:00"
  }
}

# align to genome
## Currently uses -M but GATK uses -Y and no -M
task cramToBwaMemMerge {
  input {
    File input_cram
    String base_file_name
    File ref_fasta
    File ref_fasta_index
    File ref_dict
    File? ref_alt
    File ref_amb
    File ref_ann
    File ref_bwt
    File ref_pac
    File ref_sa
    String taskDocker
    Int bwaThreads
  }
  command {
    set -eo pipefail

    java -Xms5000m -jar /usr/gitc/picard.jar \
    SamToFastq \
      INPUT=~{input_bam} \
      FASTQ=/dev/stdout \
      INTERLEAVE=true \
      NON_PF=true | \
    bwa mem \
      -K 100000000 \
      -p -v 3 -t ~{bwaThreads} -Y \
      ~{ref_fasta} /dev/stdin - 2> >(tee ~{output_bam_basename}.bwa.stderr.log >&2) | \
    gatk --java-options "-Dsamjdk.compression_level=5 -XX:-UseGCOverheadLimit -Xms12g" \
      MergeBamAlignment \
      --VALIDATION_STRINGENCY SILENT \
      --EXPECTED_ORIENTATIONS FR \
      --ATTRIBUTES_TO_RETAIN X0 \
      --ALIGNED_BAM ~{aligned_bam} \
      --UNMAPPED_BAM ~{unmapped_bam} \
      --OUTPUT ~{base_file_name}.merged.bam \
      --REFERENCE_SEQUENCE ~{ref_fasta} \
      --PAIRED_RUN true \
      --SORT_ORDER queryname \
      --IS_BISULFITE_SEQUENCE false \
      --ALIGNED_READS_ONLY false \
      --CLIP_ADAPTERS false \
      --MAX_RECORDS_IN_RAM 500000 \
      --ADD_MATE_CIGAR true \
      --MAX_INSERTIONS_OR_DELETIONS -1 \
      --PRIMARY_ALIGNMENT_STRATEGY MostDistant \
      --UNMAPPED_READ_STRATEGY COPY_TO_TAG \
      --ALIGNER_PROPER_PAIR_FLAGS true \
      --CREATE_INDEX false
  }
  output {
    File output_bam = "~{base_file_name}.aligned.bam"
  }
  runtime {
    memory: "48 GB"
    cpu: bwaThreads
    docker: taskDocker
  }
}



# Merge original input uBAM file with BWA-aligned BAM file
task MergeBamAlignment {
  input {
    File unmapped_bam
    File aligned_bam
    String base_file_name
    File ref_fasta
    File ref_fasta_index
    File ref_dict
    String taskDocker
  }
  command {
    set -eo pipefail

    gatk --java-options "-Dsamjdk.compression_level=5 -XX:-UseGCOverheadLimit -Xms12g" \
      MergeBamAlignment \
      --VALIDATION_STRINGENCY SILENT \
      --EXPECTED_ORIENTATIONS FR \
      --ATTRIBUTES_TO_RETAIN X0 \
      --ALIGNED_BAM ~{aligned_bam} \
      --UNMAPPED_BAM ~{unmapped_bam} \
      --OUTPUT ~{base_file_name}.merged.bam \
      --REFERENCE_SEQUENCE ~{ref_fasta} \
      --PAIRED_RUN true \
      --SORT_ORDER queryname \
      --IS_BISULFITE_SEQUENCE false \
      --ALIGNED_READS_ONLY false \
      --CLIP_ADAPTERS false \
      --MAX_RECORDS_IN_RAM 500000 \
      --ADD_MATE_CIGAR true \
      --MAX_INSERTIONS_OR_DELETIONS -1 \
      --PRIMARY_ALIGNMENT_STRATEGY MostDistant \
      --UNMAPPED_READ_STRATEGY COPY_TO_TAG \
      --ALIGNER_PROPER_PAIR_FLAGS true \
      --CREATE_INDEX false
  }
  output {
    File output_bam = "~{base_file_name}.merged.bam"
  }
  runtime {
    memory: "16 GB"
    cpu: 2
    docker: taskDocker
  }
}

# Generate Base Quality Score Recalibration (BQSR) model and apply it
task ApplyBaseRecalibrator {
  input {
    File input_bam
    File intervals 
    File input_bam_index
    String base_file_name
    File dbSNP_vcf
    File dbSNP_vcf_index
    Array[File] known_indels_sites_VCFs
    Array[File] known_indels_sites_indices
    File ref_dict
    File ref_fasta
    File ref_fasta_index
    String taskDocker
  }
  command {
  set -eo pipefail

  samtools index ~{input_bam}

  gatk --java-options "-Xms8g" \
      BaseRecalibrator \
      -R ~{ref_fasta} \
      -I ~{input_bam} \
      -O ~{base_file_name}.recal_data.csv \
      --known-sites ~{dbSNP_vcf} \
      --known-sites ~{sep=" --known-sites " known_indels_sites_VCFs} \
      --intervals ~{intervals} \
      --interval-padding 100 

  gatk --java-options "-Xms8g" \
      ApplyBQSR \
      -bqsr ~{base_file_name}.recal_data.csv \
      -I ~{input_bam} \
      -O ~{base_file_name}.recal.bam \
      -R ~{ref_fasta} \
      --intervals ~{intervals} \
      --interval-padding 100 

  #finds the current sort order of this bam file
  samtools view -H ~{base_file_name}.recal.bam | grep @SQ | sed 's/@SQ\tSN:\|LN://g' > ~{base_file_name}.sortOrder.txt

  }
  output {
    File recalibrated_bam = "~{base_file_name}.recal.bam"
    File recalibrated_bai = "~{base_file_name}.recal.bai"
    File sortOrder = "~{base_file_name}.sortOrder.txt"
  }
  runtime {
    memory: "36 GB"
    cpu: 2
    docker: taskDocker
  }
}



task MarkDuplicatesSpark {
  input {
    Array[File] input_bams
    String output_bam_basename
    String metrics_filename
    String taskDocker
  }
  # LAter use: --verbosity WARNING
 # Task is assuming query-sorted input so that the Secondary and Supplementary reads get marked correctly.
 # This works because the output of BWA is query-grouped and therefore, so is the output of MergeBamAlignment.
 # While query-grouped isn't actually query-sorted, it's good enough for MarkDuplicates with ASSUME_SORT_ORDER="queryname"
  command {
    gatk --java-options "-XX:+UseParallelGC -XX:ParallelGCThreads=4 -Dsamjdk.compression_level=5 -Xmx64g" \
      MarkDuplicatesSpark \
      --input ~{sep=" --input " input_bams} \
      --output ~{output_bam_basename}.bam \
      --metrics-file ~{metrics_filename} \
      --create-output-bam-index true \
      --spark-verbosity WARN \
      --verbosity WARNING \
      --optical-duplicate-pixel-distance 2500 
  }
  runtime {
    docker: taskDocker
    memory: "64 GB"
    cpu: 32
  }
  output {
    File output_bam = "~{output_bam_basename}.bam"
    File output_bai = "~{output_bam_basename}.bam.bai"
    File duplicate_metrics = "~{metrics_filename}"
  }
}

task SplitIntervals {
    input {
      File intervals
      File ref_fasta
      File ref_fasta_index
      File ref_dict
      Int scatter_count
      String taskDocker
    }
    command {
        set -eo pipefail
        mkdir interval-files1
        gatk --java-options "-Xms2g" \
          SplitIntervals \
            -R ~{ref_fasta} \
            -L ~{intervals} \
            -scatter ~{scatter_count} \
            -O interval-files1 \
            --subdivision-mode BALANCING_WITHOUT_INTERVAL_SUBDIVISION

        cp interval-files1/*.interval_list .
        ls *.interval_list > globGetAround.txt
    }
    runtime {
        docker: taskDocker
        cpu: 2
        memory: "2 GB"
    }
    output {
        Array[File] interval_files = read_lines("globGetAround.txt")
    }
}