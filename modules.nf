// Procesos llamado de varialtes germinal
process fastqc {
  cache 'lenient'
  publishDir params.out, mode:'symlink'
  
  input: 
  tuple val(sample), path(reads)
  
  output:
  path("fastqc/${sample}/*"), emit: fq_files

  script:
  """
    mkdir -p fastqc/${sample}
    fastqc -o fastqc/${sample} -f fastq -q ${reads[0]} ${reads[1]}
  """   
}

process multiqc {
  cache 'lenient'
  publishDir params.out, mode: 'symlink'
    
  input:
  val(sample_1)
  
  output:
  path("multiqc/*")   , emit: multiqc_fq_data
  
  script: 
  
  x_1 = sample_1.size()
  
  """
    echo -e "\n"
   
    mkdir -p multiqc
           
    multiqc -o multiqc/ ${params.out}/fastqc/
    
    echo -e "\nEl número de archivos encontrados en el directorio 'fastqc' es $x_1 \n"
  """  
}

process alignG {
    cache 'lenient'
    publishDir params.out + "/aligned_reads", mode:'symlink'

    input:
    tuple val(pair_id), path(reads) 
         
    output:
    tuple val(pair_id), path("${pair_id}_aligned_reads.sam"),     emit: aligned_reads_ch
    
    script:
    
    FLOWCELL="\$(zcat ${reads[0]} | head -n 1 | awk -F : '{ print \$3 }')"
    LANE="\$(zcat ${reads[0]} | head -n 1 | awk -F : '{ print \$4 }')"
    sample_id = "\$(echo ${pair_id} | cut -d'_' -f 2)"
    
    readGroup = \
	"@RG\\tID:${pair_id}\\tLB:${sample_id} '.' ${FLOWCELL} '.' ${LANE}\\tPL:${params.pl}\\tPM:${params.pm}\\tSM:${sample_id}"
	
    """     
    bwa mem \
	-K 100000000 \
	-v 3 \
	-t ${params.ncrs} \
	-Y \
	-R \"${readGroup}\" \
	${params.ref} \
	${reads[0]} \
	${reads[1]} > ${pair_id}_aligned_reads.sam
    """	
}

process mergeSam {
   cache 'lenient'
   publishDir params.out + "/merged_sam", mode:'symlink'

   input:
   tuple val(key), file(bam_files)

   output:
   tuple val(key), path("${key}_merged.sam"),     emit: merged_sam_ch

   script:
   """
   java -jar ${params.picard} MergeSamFiles \
            ${'INPUT='+bam_files.join(' INPUT=')} \\
            OUTPUT=${key}_merged.sam
   """
}

process markDuplicatesSpark {
    cache 'lenient'
    publishDir params.out + "/dedup_sorted", mode:'symlink'

    input:
    tuple val(pair_id), path(aligned_reads)

    output:
    tuple val(pair_id), path("${pair_id}_sorted_dedup.bam"),    emit: bam_for_variant_calling
    tuple val(pair_id), path("${pair_id}_dedup_metrics.txt"),   emit: dedup_qc_ch

    script:
    // --java-options "-Djava.io.tmpdir=${params.tmpdir}/${workflow.runName}/${pair_id} -Xms16G -Xmx16G "
    """   
    export JAVA_HOME="/labs/sbio/dpservs/miniconda3/envs/gatk/bin"
    export PATH=/labs/sbio/dpservs/miniconda3/envs/gatk/bin/:$PATH
    
    mkdir -p ${params.tmpdir}/${workflow.runName}/${pair_id}
    
    ${params.gatk} --java-options "-Xms16G -Xmx16G" \
	 MarkDuplicatesSpark \
	-I ${aligned_reads} \
	-M ${pair_id}_dedup_metrics.txt \
	-O ${pair_id}_sorted_dedup.bam \
        --tmp-dir ${params.tmpdir}/${workflow.runName}/${pair_id}  
    rm -r ${params.tmpdir}/${workflow.runName}/${pair_id}
    """ 
}

process getMetrics {
    cache 'lenient'
    publishDir params.out + "/metrics", mode:'symlink'

    input:
    tuple val(pair_id), path(sorted_dedup_reads)

    output:
    tuple val(pair_id), 
	  path("${pair_id}_alignment_metrics.txt"), \
	  path("${pair_id}_insert_metrics.txt"), \
	  path("${pair_id}_insert_size_histogram.pdf"), \
	  path("${pair_id}_depth_out.txt"), \
	  path("${pair_id}_cov_hist.txt"),             emit: metrics_qc_ch

    script:
	//samtools coverage -w 32 -o ${pair_id}_cov_hist.txt ${sorted_dedup_reads}
    """  
    mkdir -p metrics
    java -jar ${params.picard} \
        CollectAlignmentSummaryMetrics \
	-R ${params.ref} \
        -I ${sorted_dedup_reads} \
	-O ${pair_id}_alignment_metrics.txt
    java -jar ${params.picard} \
        CollectInsertSizeMetrics \
        -I ${sorted_dedup_reads} \
	-O ${pair_id}_insert_metrics.txt \
        -H ${pair_id}_insert_size_histogram.pdf 
    samtools depth -a ${sorted_dedup_reads} > ${pair_id}_depth_out.txt
    samtools coverage -w 32 -o ${pair_id}_cov_hist.txt ${sorted_dedup_reads}
    """
}

process haplotypeCaller {
    cache 'lenient'
    
    input:
    tuple val(pair_id), path(input_bam)

    output:
    tuple val(pair_id), path("${pair_id}_raw_variants.vcf"), emit: hc_output

    script:
    """
    ${params.gatk} HaplotypeCaller \
	-R ${params.ref} \
	-I $input_bam \
	-O ${pair_id}_raw_variants.vcf \
    """
}

process selectVariants {
    cache 'lenient'
    
    input:
    tuple val(pair_id), path(raw_variants)

    output:
    tuple val(pair_id), path("${pair_id}_raw_snps_0.vcf"),    emit: raw_snps_ch
    tuple val(pair_id), path("${pair_id}_raw_indels_0.vcf"),  emit: raw_indels_ch

    script:
    """
    ${params.gatk} SelectVariants \
	-R ${params.ref} \
	-V ${raw_variants} \
	-select-type SNP \
	-O ${pair_id}_raw_snps_0.vcf
    ${params.gatk} SelectVariants \
        -R ${params.ref} \
        -V ${raw_variants} \
        -select-type INDEL \
        -O ${pair_id}_raw_indels_0.vcf
    """
}

process filterSnps {
    cache 'lenient'
    publishDir params.out + "/bootstrapping", mode:'symlink'
    
    input:
    tuple val(pair_id), path(raw_snps)

    output:
    tuple val(pair_id), \
    path("${pair_id}_filtered_snps.vcf"), \
    path("${pair_id}_filtered_snps.vcf.idx"),     emit: filtered_snps
    
    script:
    """
    ${params.gatk} VariantFiltration \
	-R ${params.ref} \
	-V ${raw_snps} \
	-O ${pair_id}_filtered_snps.vcf \
	-filter-name "QD_filter" -filter "QD < 2.0" \
	-filter-name "FS_filter" -filter "FS > 60.0" \
	-filter-name "MQ_filter" -filter "MQ < 40.0" \
	-filter-name "SOR_filter" -filter "SOR > 4.0" \
	-filter-name "MQRankSum_filter" -filter "MQRankSum < -12.5" \
	-filter-name "ReadPosRankSum_filter" -filter "ReadPosRankSum < -8.0"
    """
}

process filterIndels {
    cache 'lenient'
    publishDir params.out + "/bootstrapping", mode:'symlink'
    
    input:
    tuple val(pair_id), path(raw_indels) 

    output:
    tuple val(pair_id), \
	  path("${pair_id}_filtered_indels.vcf"), \
	  path("${pair_id}_filtered_indels.vcf.idx"), emit: filtered_indels

    script:
    """
    ${params.gatk} VariantFiltration \
        -R ${params.ref} \
        -V ${raw_indels} \
        -O ${pair_id}_filtered_indels.vcf \
	-filter-name "QD_filter" -filter "QD < 2.0" \
	-filter-name "FS_filter" -filter "FS > 200.0" \
	-filter-name "SOR_filter" -filter "SOR > 10.0"
    """
}

process bqsr {
    cache 'lenient'
    publishDir params.out + "/bqsr", mode:'copy'

    input:
    // tuple val(pair_id), path(input_bam)
    // tuple val(sample_snps), path(filtered_snps), path(filtered_snps_index)
    // tuple val(sample_indels), path(filtered_indels), path(filtered_indels_index)
    tuple val(sample), path(input_bam), path(filtered_snps), path(filtered_snps_index), path(filtered_indels), path(filtered_indels_index)  

    output:    
    tuple val(pair_id), path("bqsr/${pair_id}_recal_data.table"), path("bqsr/${pair_id}_post_recal_data.table"), emit: analyze_covariates
    tuple val(pair_id), path("bqsr/${pair_id}_recal.bam"),     emit: recalibrated_bam
    path("bqsr/${pair_id}_recal.bai")

    script:
    """
    ${params.gatk} SelectVariants \
	--exclude-filtered \
	-V ${filtered_snps} \
	-O ${sample}_bqsr_snps.vcf
    ${params.gatk} SelectVariants \
        --exclude-filtered \
        -V ${filtered_indels} \
        -O ${sample}_bqsr_indels.vcf
        
    export JAVA_HOME="/labs/sbio/dpservs/miniconda3/envs/gatk/bin"
    export PATH=/labs/sbio/dpservs/miniconda3/envs/gatk/bin/:$PATH
    
    ${params.gatk} BaseRecalibrator \
	-R ${params.ref} \
	-I ${input_bam} \
	--known-sites ${sample}_bqsr_snps.vcf \
	--known-sites ${sample}_bqsr_indels.vcf \
	-O ${sample}_recal_data.table
    ${params.gatk} ApplyBQSR \
        -R ${params.ref} \
        -I ${input_bam} \
        -bqsr ${sample}_recal_data.table \
        -O ${sample}_recal.bam
    ${params.gatk} BaseRecalibrator \
        -R ${params.ref} \
	-I ${sample}_recal.bam \
        --known-sites ${sample}_bqsr_snps.vcf \
	--known-sites ${sample}_bqsr_indels.vcf \
	-O ${sample}_post_recal_data.table
    """	
}

process recal_bqsr {
    cache 'lenient'
    publishDir params.out + "/bqsr", mode:'copy'

    input:
    tuple val(sample), path(reads)

    output:
    tuple val(sample), path("${sample}_recalibration_data.table"), path("${sample}_post_recal_data.table"), emit: analyze_covariates
    tuple val(sample), path("${sample}_recalibrated.bam"),     emit: recalibrated_bam
    path("${sample}_recalibrated.bai")

    script:
    """
    export JAVA_HOME="/labs/sbio/dpservs/miniconda3/envs/gatk/bin"
    export PATH=/labs/sbio/dpservs/miniconda3/envs/gatk/bin/:$PATH
    
    ${params.gatk} BaseRecalibrator \
    -I ${reads} \
    -R ${params.ref} \
    --known-sites ${params.refdir}/1000G_phase1.snps.high_confidence.hg38.vcf.gz \
    --known-sites ${params.refdir}/Homo_sapiens_assembly38.dbsnp138.vcf \
    --known-sites ${params.refdir}/Homo_sapiens_assembly38.known_indels.vcf.gz \
    --known-sites ${params.refdir}/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz \
    -O ${sample}_recalibration_data.table \
    --tmp-dir ${params.tmpdir} 
    
    ${params.gatk} ApplyBQSR \
    -R ${params.ref} \
    -I ${reads} \
    -bqsr ${sample}_recalibration_data.table \
    -O ${sample}_recalibrated.bam \
    --tmp-dir ${params.tmpdir}

    ${params.gatk} BaseRecalibrator \
    -R ${params.ref} \
    -I ${sample}_recalibrated.bam \
    --known-sites ${params.refdir}/1000G_phase1.snps.high_confidence.hg38.vcf.gz \
    --known-sites ${params.refdir}/Homo_sapiens_assembly38.dbsnp138.vcf \
    --known-sites ${params.refdir}/Homo_sapiens_assembly38.known_indels.vcf.gz \
    --known-sites ${params.refdir}/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz \
    -O ${sample}_post_recal_data.table \
    --tmp-dir ${params.tmpdir}
    """ 
}

process analyzeCovariates{
    cache 'lenient'
    publishDir params.out + "/bqsr", mode:'symlink'

    input:
    tuple val(pair_id), path(recal_table), path(post_recal_table)

    output:
    tuple val(pair_id), path("${pair_id}_recalibration_plots.pdf"), emit: analyzed_covariates_ch

    script:
    """
    export JAVA_HOME="/labs/sbio/dpservs/miniconda3/envs/gatk/bin"
    export PATH=/labs/sbio/dpservs/miniconda3/envs/gatk/bin/:$PATH
    
    ${params.gatk} AnalyzeCovariates \
	-before ${recal_table} \
	-after ${post_recal_table} \
	-plots ${pair_id}_recalibration_plots.pdf
     """
}

process haplotypeCallerERC {
    cache 'lenient'
    publishDir params.out + "/RAW_gvcfs", mode:'symlink'
    
    input:
    tuple val(pair_id), path(input_bam)

    output:
    tuple val(pair_id), path("${pair_id}_raw_variants.g.vcf.gz"), emit: hc_erc_out
    path("${pair_id}_raw_variants.g.vcf.*")

    script:
    """
    ${params.gatk} HaplotypeCaller \
	-R ${params.ref} \
	-I $input_bam \
	-O ${pair_id}_raw_variants.g.vcf.gz \
	-max-mnp-distance 0 \
	-ERC GVCF 
    """
}

process genomicsDBimport {
    cache 'lenient'
    publishDir params.out + "/Genomicsdb", mode:'symlink'
    
    input:
    val(sample_map)
    val(project_id)

    output:
    tuple val(project_id), path("${project_id}_database"), emit: genomics_db

    script:
    """
    echo -e $sample_map | sed 's/^.//' | sed 's/.\$//' | sed 's/],/-/g' | tr '-' "\n" | tr ',' "	" | tr -d "] [" >> cohort.sample_map
    
    mkdir -p ${params.tmpdir}/${workflow.runName}/temp_genomicsdb
    
    ${params.gatk} --java-options "-Xmx16g -Xms16g" GenomicsDBImport \
       --genomicsdb-workspace-path ${project_id}_database \
       --batch-size ${params.batchsize} \
       -L ${params.interval_list} \
       --sample-name-map cohort.sample_map \
       --tmp-dir ${params.tmpdir}/${workflow.runName}/temp_genomicsdb \
       --reader-threads ${params.GI_cores}
       
    rm -r ${params.tmpdir}/${workflow.runName}/temp_genomicsdb
    """
}

process genotypeGVCFs {
    cache 'lenient'
    publishDir params.out + "/join_vcfs", mode:'symlink'
    
    input:
    tuple val(project_id), path(database)

    output:
    tuple val(project_id), path("${project_id}_raw_variants.vcf.gz"), emit: gvcfs_out
    path("${project_id}_raw_variants.vcf.gz.tbi")

    script:
    """
    ${params.gatk} --java-options "-Xmx8g" GenotypeGVCFs \
   -R ${params.ref} \
   -V gendb://${database} \
   -O ${project_id}_raw_variants.vcf.gz
    """
}

process selectVariantsConjunto {
    cache 'lenient'
    publishDir params.out + "/join_vcfs", mode:'symlink'
    
    input:
    tuple val(project_id), val(variants)

    output:
    tuple val(project_id), path("${project_id}_raw_snps.vcf"),    emit: snps_ch
    tuple val(project_id), path("${project_id}_raw_indels.vcf"),  emit: indels_ch

    script:
    """
    ${params.gatk} SelectVariants \
	-R ${params.ref} \
	-V ${variants} \
	-select-type SNP \
	-O ${project_id}_raw_snps.vcf
    ${params.gatk} SelectVariants \
        -R ${params.ref} \
        -V ${variants} \
        -select-type INDEL \
        -O ${project_id}_raw_indels.vcf
    """
}

process VQSRsnps {
    cache 'lenient'
    publishDir params.out + "/join_vcfs", mode:'symlink'
    
    input:
    tuple val(project_id), path(raw_snps)

    output:
    tuple val(project_id), path("${project_id}_filtered_snps.vcf"), path("${project_id}_filtered_snps.vcf.idx"),    emit: snps_filt_ch

    script:
    """   
    ${params.gatk} VariantRecalibrator \
   -R ${params.ref} \
   -V ${raw_snps} \
   --resource:hapmap,known=false,training=true,truth=true,prior=15.0 ${params.refdir}/hapmap_3.3.hg38.vcf.gz \
   --resource:omni,known=false,training=true,truth=false,prior=12.0 ${params.refdir}/1000G_omni2.5.hg38.vcf.gz \
   --resource:1000G,known=false,training=true,truth=false,prior=10.0 ${params.refdir}/1000G_phase1.snps.high_confidence.hg38.vcf.gz \
   --resource:dbsnp,known=true,training=false,truth=false,prior=2.0 ${params.refdir}/Homo_sapiens_assembly38.dbsnp138.vcf \
   -an QD -an MQ -an MQRankSum -an ReadPosRankSum -an FS -an SOR \
   -mode SNP \
   -O ${project_id}_snps_output.recal \
   --tranches-file ${project_id}_snps_output.tranches \
   --rscript-file ${project_id}_snps_output.plots.R
   
   ${params.gatk} ApplyVQSR \
   -R ${params.ref} \
   -V ${raw_snps} \
   -O ${project_id}_filtered_snps.vcf \
   -ts-filter-level 99.5 \
   --tranches-file ${project_id}_snps_output.tranches \
   --recal-file ${project_id}_snps_output.recal \
   -mode SNP \
   --create-output-variant-index true  
    """
}

process VQSRindels {
    cache 'lenient'
    publishDir params.out + "/join_vcfs", mode:'symlink'
    
    input:
    tuple val(project_id), path(raw_indels)

    output:
    tuple val(project_id), path("${project_id}_filtered_indels.vcf"), path("${project_id}_filtered_indels.vcf.idx"),   emit: indels_filt_ch

    script:
    """   
    ${params.gatk} VariantRecalibrator \
   -R ${params.ref} \
   -V ${raw_indels} \
   -resource:mills,known=false,training=true,truth=true,prior=12.0 ${params.refdir}/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz \
   -resource:dbsnp,known=true,training=false,truth=false,prior=2.0 ${params.refdir}/Homo_sapiens_assembly38.dbsnp138.vcf \
   -an QD -an MQRankSum -an ReadPosRankSum -an FS -an SOR -an DP \
   -mode INDEL \
   -O ${project_id}_indels_output.recal \
   --tranches-file ${project_id}_indels_output.tranches \
   --rscript-file ${project_id}_indels_output.plots.R
   
   ${params.gatk} ApplyVQSR \
   -R ${params.ref} \
   -V ${raw_indels} \
   -O ${project_id}_filtered_indels.vcf \
   -ts-filter-level 99.0 \
   --tranches-file ${project_id}_indels_output.tranches \
   --recal-file ${project_id}_indels_output.recal \
   -mode INDEL \
   --create-output-variant-index true
    """
}

process snpEff {
    cache 'lenient'
    publishDir params.out + "/snpeff_vcf", mode:'symlink'    

    input:
    tuple val(project_id), path(filtered_snps)  , path(filtered_snps_index)
    tuple val(project_id), path(filtered_indels), path(filtered_indels_index) 

    output:
    tuple val("${project_id}_snpEff"), path("${project_id}_filtered_snps.ann.vcf"), path("${project_id}_filtered_indels.ann.vcf") , emit: snpeff_ch

    script:
    """
    java -jar ${params.snpEff} -v \
	-dataDir ${params.snpeff_data} \
	${params.snpeff_db} \
	${filtered_snps} > ${project_id}_filtered_snps.ann.vcf

    java -jar ${params.snpEff} -v \
        -dataDir ${params.snpeff_data} \
        ${params.snpeff_db} \
        ${filtered_indels} > ${project_id}_filtered_indels.ann.vcf
    """
}

process annovarG {
    cache 'lenient'
    publishDir params.out, mode:'symlink' 

    input:
    tuple val(project_id), path(filtered_snps)  , path(filtered_snps_index)
    tuple val(project_id), path(filtered_indels), path(filtered_indeles_index)

    output:
    tuple val("${project_id}_annovar"), path("annovar/${project_id}_annovar_snps*.vcf"), path("annovar/${project_id}_annovar_indels*.vcf"), emit: annovar_ch
    path("annovar/*")
    
    script:
    """
    mkdir -p annovar
    perl ${params.annovar}/table_annovar.pl ${filtered_snps} ${params.annovar}/humandb/ --buildver hg38 -out annovar/${project_id}_annovar_snps \
    -remove -protocol refGene,ensGene,avsnp150,clinvar_20220320,gnomad_genome,cosmic70,dbnsfp33a -operation g,g,f,f,f,f,f -nastring . -vcfinput

    perl ${params.annovar}/table_annovar.pl ${filtered_indels} ${params.annovar}/humandb/ --buildver hg38 -out annovar/${project_id}_annovar_indels \
    -remove -protocol refGene,ensGene,avsnp150,clinvar_20220320,gnomad_genome,cosmic70,dbnsfp33a -operation g,g,f,f,f,f,f -nastring . -vcfinput
    """
}

process splitVCFs {
    cache 'lenient'
    publishDir params.out + "/vcfs_persample" , mode:'copy'

    input:
    tuple val(project_id), path(project_vcf_snps), path(project_vcf_indels)

    output:
    path("*.ann.vcf"),    emit: ann_vcf_persample

    script:
    """
    for sample in `bcftools query -l ${project_vcf_snps}`
    do
    ${params.gatk} SelectVariants -R ${params.ref} -V ${project_vcf_snps} -O \${sample}_filtered_snps_${project_id}.ann.vcf -sn \${sample} --exclude-non-variants true
    done

    for sample in `bcftools query -l ${project_vcf_indels}`
    do
    ${params.gatk} SelectVariants -R ${params.ref} -V ${project_vcf_indels} -O \${sample}_filtered_indels_${project_id}.ann.vcf -sn \${sample} --exclude-non-variants true
    done
    """
}

// Procesos llamado de variantes somático
process fastq_to_sam {
    cache 'lenient'
    publishDir params.out + "/uBams", mode:'symlink'

    input:
    tuple val(sample), val(RG), val(PU), path(read1), path(read2)
         
    output:
    tuple val(sample), path("${sample}_fastqtosam.bam"),  emit: fastq_to_sam_ch
	
    script:   
    """
    export JAVA_HOME="/labs/sbio/dpservs/miniconda3/envs/gatk/bin"
    export PATH=/labs/sbio/dpservs/miniconda3/envs/gatk/bin/:$PATH
    
    mkdir -p uBams

   java -Xmx8G -jar ${params.picard} FastqToSam \
    -FASTQ ${read1} \
    -FASTQ2 ${read2} \
    -OUTPUT ${sample}_fastqtosam.bam \
    -READ_GROUP_NAME ${RG} \
    -SAMPLE_NAME ${sample} \
    -LIBRARY_NAME ${sample} \
    -PLATFORM_UNIT ${PU} \
    -PLATFORM illumina 
    """
}

process mark_duplicates {
    cache 'lenient'
    publishDir params.out + "/MarkedDuplicates", mode:'symlink'

    input:
    tuple val(sample), path(reads) 
         
    output:
    tuple val(sample), path("${sample}_markilluminaadapters.bam"),          emit: mark_duplicates_bam
    tuple val(sample), path("${sample}_markilluminaadapters_metrics.txt"),  emit: mark_duplicates_metrics

    script:
    """
    export JAVA_HOME="/labs/sbio/dpservs/miniconda3/envs/gatk/bin"
    export PATH=/labs/sbio/dpservs/miniconda3/envs/gatk/bin/:$PATH
    
    java -Xmx8G -jar ${params.picard} MarkIlluminaAdapters \
      -I ${reads} \
      -O ${sample}_markilluminaadapters.bam \
      -M ${sample}_markilluminaadapters_metrics.txt \
      --TMP_DIR ${params.tmp}
    """
}

process sam_to_fastq {
    cache 'lenient'
    publishDir params.out + "/sam_to_fq", mode:'symlink'

    input:
    tuple val(sample), path(reads) 
         
    output:
    tuple val(sample), path("${sample}_samtofastq.fastq"),     emit: sam_to_fastq_ch
	
    script:
    """
    export JAVA_HOME="/labs/sbio/dpservs/miniconda3/envs/gatk/bin"
    export PATH=/labs/sbio/dpservs/miniconda3/envs/gatk/bin/:$PATH
    
    mkdir -p Alignments
    java -Xmx16G -jar ${params.picard} SamToFastq \
      -I ${reads} \
	    -FASTQ ${sample}_samtofastq.fastq \
	    -CLIPPING_ATTRIBUTE XT \
      -CLIPPING_ACTION 2 \
      -INTERLEAVE true \
      -NON_PF true \
      --TMP_DIR ${params.tmp}
    """
}

process alignS {
    cache 'lenient'
    publishDir params.out +"/Alignments", mode:'symlink'

    input:
    tuple val(sample), path(reads) 
         
    output:
    tuple val(sample), path("${sample}_bwa_mem.sam"),     emit: aligned_reads
	
    script:
    """   
    bwa mem -M -t 4 -p ${params.ref} ${reads} >  ${sample}_bwa_mem.sam
    """
}

process merge_bam_alignment {
  cache 'lenient'
  publishDir params.out + "/merge_bam_algn", mode:'copy'
  
  input: 
  tuple val(sample), path(reads_1), path(reads_2)
  
  output:
  tuple val(sample), path("${sample}_merged.bam"),     emit: bam_alignment_ch

  script:
  """
    export JAVA_HOME="/labs/sbio/dpservs/miniconda3/envs/gatk/bin"
    export PATH=/labs/sbio/dpservs/miniconda3/envs/gatk/bin/:$PATH
    
  java -Xmx16G -jar ${params.picard} MergeBamAlignment \
      -ALIGNED ${reads_1} \
      -UNMAPPED ${reads_2} \
      -O ${sample}_merged.bam \
      -R ${params.ref} \
      --CREATE_INDEX true \
      --ADD_MATE_CIGAR true \
      --CLIP_ADAPTERS false \
      --CLIP_OVERLAPPING_READS true \
      --INCLUDE_SECONDARY_ALIGNMENTS true \
      --MAX_INSERTIONS_OR_DELETIONS -1 \
      --PRIMARY_ALIGNMENT_STRATEGY MostDistant \
      --ATTRIBUTES_TO_RETAIN XS \
      --TMP_DIR ${params.tmp}
  """   
}

process mutect2forPanelofNormals {
    cache 'lenient'
    publishDir params.out + "/vcfsforPON", mode:'symlink'
    
    input:
    tuple val(pair_id), path(input_bam)

    output:
    tuple val(pair_id), path("${pair_id}_for_pon.vcf.gz"), emit: mtf_PON_out

    script:
    // --genotype-germline-sites true
    """
   ${params.gatk} Mutect2 \
   -R ${params.ref} \
   -I ${input_bam} \
   -max-mnp-distance 0 \
   -tumor ${pair_id} \
   --germline-resource ${params.onlygnomad} \
   -O ${pair_id}_for_pon.vcf.gz   
    """
}

process createSomaticPanelofNormals {
    cache 'lenient'
    publishDir params.out + "/panelofNormals", mode:'copy'
       
    input:
    val(input_db)

    output:
    path("m_pon.vcf.gz"), emit: PON_out
    path("m_pon.vcf.gz.tbi")

    script:
    """
    ${params.gatk} CreateSomaticPanelOfNormals -R ${params.ref} --germline-resource ${params.onlygnomad} --output m_pon.vcf.gz -V gendb://${input_db}
    """
}

process mutect2 {
    cache 'lenient'
    publishDir params.out + "/unfilteredVCFs", mode:'symlink'

    input:
    tuple val(pair_id), path(input_bam)

    output:
    tuple val(pair_id), path("${pair_id}_unfiltered.vcf"), path("${pair_id}_read-orientation-model.tar.gz"), emit: unfilt_vcf
    path("${pair_id}_unfiltered.vcf.stats")

    script:
    """
   ${params.gatk} Mutect2 \
   -R ${params.ref} \
   -L ${params.interval_list} \
   -I ${input_bam} \
   --germline-resource ${params.onlygnomad} \
   -pon ${params.panelnormales} \
   --genotype-germline-sites true \
   --genotype-pon-sites true \
   --f1r2-tar-gz ${pair_id}_f1r2.tar.gz \
   -O ${pair_id}_unfiltered.vcf

   ${params.gatk} LearnReadOrientationModel -I ${pair_id}_f1r2.tar.gz -O ${pair_id}_read-orientation-model.tar.gz
    """
}

process calculateContamination {
    cache 'lenient'
    publishDir params.out + "/contamination_tables", mode:'symlink'

    input:
    tuple val(tumor_id), path(input_bam)

    output:
    tuple val(tumor_id), path("${tumor_id}_segments.tsv"), path("${tumor_id}_calculatecontamination.table"),  emit: cont_tables

    script:
    """
    ${params.gatk} GetPileupSummaries \
    -I ${input_bam} \
    -V ${params.common_biallelic}  \
    -L ${params.common_biallelic}  \
    -O ${tumor_id}_getpileupsummaries.table

   ${params.gatk} CalculateContamination \
        -I ${tumor_id}_getpileupsummaries.table \
        -tumor-segmentation ${tumor_id}_segments.tsv \
        -O ${tumor_id}_calculatecontamination.table
    """
}

process filterMutectCalls {
    cache 'lenient'
    publishDir params.out + "/filteredVCF" , mode:'copy'

    input:
    tuple val(tumor_id), val(input_vcf), path(orient_model), path(seg_table), path(cont_table) 

    output:
    tuple val(tumor_id), path("${tumor_id}_filtered.vcf"),    emit: filt_vcf
    path("${tumor_id}_filtered.vcf.*")

    script:
    """
    ${params.gatk} FilterMutectCalls \
        -R ${params.ref} \
        -V ${input_vcf} \
        -O ${tumor_id}_filtered.vcf \
        --tumor-segmentation ${seg_table} \
        --contamination-table ${cont_table} \
        --ob-priors ${orient_model}
    """
}

process filterSVCF {
    cache 'lenient'
    publishDir params.out + "/filteredVCF" , mode:'symlink'

    input:
    tuple val(tumor_id), val(input_vcf)

    output:
    tuple val(tumor_id), path("${tumor_id}_filtered_dp10_af10.vcf"),    emit: filt_vcf

    script:
    """
    grep -v "^##" ${input_vcf} | grep -E 'QUAL|PASS' | awk '{if (NR==1) {print \$0;} {split(\$8,a,";"); if (a[3]~"DP") \
    {split(a[3],b,"="); if (b[2]>10) {split(\$10,c,":"); split(c[2],d,","); if (d[2]>10) {print \$0}}}}}' > ${tumor_id}_filtered_dp10_af10.vcf
    """
}

process annovarS {
    cache 'lenient'
    publishDir params.out, mode:'symlink' 
    
    input:
    tuple val(pair_id), path(filtered_snps)

    output:
    path("annovar_snp/*"), emit: annovar_out
    
    script:
    """
    mkdir -p annovar_snp
    perl ${params.annovar}/table_annovar.pl ${filtered_snps} ${params.annovar}/humandb/ --buildver hg38 -out annovar_snp/${pair_id}_annovar \
    -remove -protocol refGene,ensGene,avsnp150,clinvar_20220320,gnomad_genome,cosmic70,dbnsfp33a -operation g,g,f,f,f,f,f -nastring . -vcfinput
    """
}

// Procesos cuantificación y análisis de expresión diferencial 
process multiqc_DEA {
  cache 'lenient'
  publishDir params.out, mode: 'symlink'
  debug true
    
  input:
  val(sample_1)
  val(sample_2)
  
  output:
  path("multiqc/fastq/*")   , emit: multiqc_fq_data
  path("multiqc/trim_g/*")  , emit: multiqc_tg_data
  
  script: 
  
  x_1 = sample_1.size()
  x_2 = sample_2.size() * 4
  
  """
    echo -e "\n"
   
    mkdir -p multiqc
           
    multiqc -o multiqc/fastq ${params.out}/fastqc/
    multiqc -o multiqc/trim_g ${params.out}/trimming_files/
    
    echo -e "\nEl número de archivos en el directorio 'fastqc' es $x_1 \n"
    echo -e "El número de archivos en el directorio 'trimming_files' es $x_2 \n"
  """  
}

process trim_Galore {
  cache 'lenient'
  publishDir params.out, mode: 'symlink'

  input:
  tuple val(sample), path(reads)

  output:
  tuple val(sample), path("trimming_files/tg_${sample}/*_trimmed.fq.gz")      , emit: trim_fq
  tuple val(sample), path("trimming_files/tg_${sample}/*report.txt")          , emit: trim_report
  tuple val(sample), path("trimming_files/tg_${sample}/*_trimmed_fastqc.html"), emit: trim_html
  tuple val(sample), path("trimming_files/tg_${sample}/*_trimmed_fastqc.zip") , emit: trim_zip
  
  script:
  """
    mkdir -p trimming_files
    trim_galore -o trimming_files/tg_${sample} --fastqc ${reads[0]} ${reads[1]}   
  """
}

process kallisto {
  cache 'lenient'
  publishDir params.out, mode: 'symlink'

  input:
  tuple val(sample), val(reads)
  
  output:
  tuple val(sample), path("kallisto_quants/${sample}/abundance.h5") , emit: abundance_h5
  tuple val(sample), path("kallisto_quants/${sample}/abundance.tsv"), emit: abundance_tsv
  tuple val(sample), path("kallisto_quants/${sample}/run_info.json"), emit: run_info  


  script:
  """
    mkdir -p kallisto_quants
    kallisto quant -i ${params.k_index} -o kallisto_quants/${sample} --rf-stranded ${reads[0]} ${reads[1]}
  """
}

process tximport_q {
  cache 'lenient'
  publishDir params.out, mode: 'symlink'
  debug true
  
  input:
  val(sample_k)

  output:
  path("results/${params.mcounts}")            , emit: mcounts
  path("results/${params.mcounts_tpm}")        , emit: mcounts_tpm
  path("*.log")                                , emit: R_sesion_info

  script:

  y_1 = sample_k.size() / 2
   
  """
   echo -e "\nEl número de archivos a procesar es: $y_1\n"
   
   mkdir -p results
  
   Rscript ${params.rscript_q_dir} \
   --working_dir ${params.out} \
   --sample_info ${params.csv_info} \
   --dir_quants "kallisto_quants" \
   --gtf_file ${params.gtf_file} \
   --countsmat results/${params.mcounts} \
   --countpm results/${params.mcounts_tpm}
   
   echo -e "\nCUANTIFICACIÓN COMPLETADA\n" 
  """
}

process tximport_deseq2 {
  cache 'lenient'
  publishDir params.out, mode: 'symlink'
  debug true
  
  input:
  val(sample_k)

  output:
  path("results/${params.pca_plot_name}")      , emit: pca_plot
  path("results/${params.heatmap_name}")       , emit: heatmap_f
  path("results/${params.volcano_plot_name}")  , emit: volcano_plot
  path("results/${params.results_c_file}")     , emit: results_c
  path("results/${params.results_file_name}")  , emit: results
  path("results/${params.resultsf_file_name}") , emit: results_f
  path("results/${params.mcounts}")            , emit: mcounts
  path("results/${params.mcounts_tpm}")        , emit: mcounts_tpm
  //path("results/${params.list_gt}")           , emit: gtl_n
  path("*.log")                                , emit: R_sesion_info

  script:

  y_1 = sample_k.size() / 2
   
  """
   echo -e "\nEl número de archivos a procesar es: $y_1\n"
   
   mkdir -p results
  
   Rscript ${params.rscript_DEA_dir} \
           --working_dir ${params.out} \
           --sample_info ${params.csv_info} \
           --dir_quants "kallisto_quants" \
           --gtf_file ${params.gtf_file} \
           --condition1 ${params.condition_1} \
           --condition2 ${params.condition_2} \
           --Log2FC_th ${params.th_l2fc} \
           --p_adj_th ${params.th_padj} \
           --outdir_pca results/${params.pca_plot_name} \
           --out_p_hm results/${params.heatmap_name} \
           --outdir_vp results/${params.volcano_plot_name} \
           --outdir_cvs results/${params.results_c_file} \
           --outres_cvs results/${params.results_file_name} \
           --outfilt_cvs results/${params.resultsf_file_name} \
           --countsmat results/${params.mcounts} \
           --countpm results/${params.mcounts_tpm} \
           --data_b ${params.data_set} \
           --onto ${params.onto_type} \
           --gtl_o results/${params.list_gt}
   
   echo -e "\nANÁLISIS COMPLETO\n" 
  """
}

// Procesos para el pipeline de pureCN
process coverageNormal {
    cache 'lenient'
    publishDir params.out + "/coverage/normal" , mode:'copy'
  
    input:
    tuple val(sample), val(input_bam)

    output:
    path("${sample}/*_coverage.txt.gz")       ,    emit: conv
    path("${sample}/*_coverage_loess.txt.gz") ,    emit: conv_loss
    path("${sample}/*")

    script:
    """
    mkdir ${sample}
    
    Rscript ${params.purecn}/Coverage.R \
    --out-dir ${sample} \
    --bam ${input_bam} \
    --intervals ${params.int_listN}
    """
}

process coverageTumor {
    cache 'lenient'
    publishDir params.out + "/coverage/tumor" , mode:'copy'

    input:
    tuple val(sample), val(input_bam)

    output:
    tuple val(sample), path("${sample}/*_coverage.txt.gz")       ,    emit: conv
    tuple val(sample), path("${sample}/*_coverage_loess.txt.gz") ,    emit: conv_loss
    path("${sample}/*")

    script:
    """
    mkdir ${sample}

    Rscript ${params.purecn}/Coverage.R \
    --out-dir ${sample} \
    --bam ${input_bam} \
    --intervals ${params.int_listT}
    """
}

process normalDB {
    cache 'lenient'
    publishDir params.out + "/normalDB" , mode:'copy'

    input:
    val(path_to_files)

    output:
    tuple path("normalDB_files/normalDB_hg38.rds"), path("normalDB_files/low_coverage_targets_hg38.bed"), emit: normal_files
    path("normalDB_files/*")

    script:
    //--normal-panel ${params.pondb} \
    """
    echo -e $path_to_files | sed 's/^.//' | sed 's/.\$//' | tr ',' "\n" | tr -d "] [" >> normal_coverages.list

    mkdir normalDB_files

    Rscript ${params.purecn}/NormalDB.R \
    --out-dir normalDB_files \
    --coverage-files normal_coverages.list \
    --genome hg38
    """
}

process pureCN {
    cache 'lenient'
    publishDir params.out + "/pureCN" , mode:'copy'

    input:
    tuple val(sample), path(vcf_file), path(stats_file), path(coverage_loess)
    tuple path(normalDB), path(mapping_bias)

    output:
    path("${sample}/*")

    script:
    // --mapping-bias-file ${mapping_bias} \
    """   
    mkdir ${sample}
    
    Rscript ${params.purecn}/PureCN.R \
    --out ${sample} \
    --tumor ${coverage_loess} \
    --sampleid ${sample} \
    --vcf ${vcf_file} \
    --stats-file ${stats_file} \
    --fun-segmentation PSCBS \
    --normaldb ${normalDB} \
    --intervals ${params.int_listT} \
    --genome hg38 \
    --model betabin \
    --force --post-optimize --seed 123
    """
}

// Modulos llamado de variantes RNS-seq (vcr)

process star {
  cache 'lenient'
  publishDir params.out +"/alignments", mode: 'copy'
   
  input:
    tuple val(sample), path(reads) 

  output:
    tuple val(sample), path("${sample}_*.sam"), emit: aligned_reads_ch
    //path("/*")

  script:
    FLOWCELL="\$(zcat ${reads[0]} | head -n 1 | awk -F : '{ print \$3 }')"
    LANE="\$(zcat ${reads[0]} | head -n 1 | awk -F : '{ print \$4 }')"
    sample_id ="\$(echo ${sample} | cut -d'_' -f 2)"
    
    readGroup = \
	"@RG\\tID:${sample}\\tLB:${sample_id} '.' ${FLOWCELL} '.' ${LANE}\\tPL:${params.pl}\\tPM:${params.pm}\\tSM:${sample_id}"
  """
    STAR --runMode alignReads \
       --genomeDir ${params.ref_dir} \
       --runThreadN ${params.cores} \
       --readFilesIn ${reads[0]} ${reads[1]} \
       --outFileNamePrefix ${sample}"_" \
       --outReadsUnmapped None \
       --twopassMode Basic \
       --twopass1readsN -1 \
       --readFilesCommand "zcat" \
       --outSAMattrRGline ID: \"${readGroup}\" 
  """
}

process haplotypeCaller_vcr {
    cache 'lenient'
    
    input:
    tuple val(pair_id), path(input_bam)
    val(round)

    output:
    tuple val(pair_id), path("${pair_id}_raw_variants_${round}.vcf"), emit: hc_output

    script:
    """
    ${params.gatk} HaplotypeCaller \
	-R ${params.ref} \
	-I ${input_bam} \
	-O ${pair_id}_raw_variants_${round}.vcf \
    """
}

process selectVariants_vcr {
    cache 'lenient'
    
    input:
    tuple val(pair_id), path(raw_variants)
    val(round)

    output:
    tuple val(pair_id), path("${pair_id}_raw_snps_${round}.vcf"),    emit: raw_snps_ch
    tuple val(pair_id), path("${pair_id}_raw_indels_${round}.vcf"),  emit: raw_indels_ch

    script:
    """
    ${params.gatk} SelectVariants \
	-R ${params.ref} \
	-V ${raw_variants} \
	-select-type SNP \
	-O ${pair_id}_raw_snps_${round}.vcf
    ${params.gatk} SelectVariants \
        -R ${params.ref} \
        -V ${raw_variants} \
        -select-type INDEL \
        -O ${pair_id}_raw_indels_${round}.vcf
    """
}

process filterSnps_vcr {
    cache 'lenient'
    publishDir params.out + "/filtered_snps" , mode:'symlink'
    
    input:
    tuple val(pair_id), path(raw_snps)
    val(round)

    output:
    tuple val(pair_id), \
    path("${pair_id}_filtered_snps_${round}.vcf"), \
    path("${pair_id}_filtered_snps_${round}.vcf.idx"),     emit: filtered_snps
    
    script:
    """
    ${params.gatk} VariantFiltration \
	-R ${params.ref} \
	-V ${raw_snps} \
	-O ${pair_id}_filtered_snps_${round}.vcf \
	-filter-name "QD_filter" -filter "QD < 2.0" \
	-filter-name "FS_filter" -filter "FS > 60.0" \
	-filter-name "MQ_filter" -filter "MQ < 40.0" \
	-filter-name "SOR_filter" -filter "SOR > 4.0" \
	-filter-name "MQRankSum_filter" -filter "MQRankSum < -12.5" \
	-filter-name "ReadPosRankSum_filter" -filter "ReadPosRankSum < -8.0"
    """
}

process filterIndels_vcr {
    cache 'lenient'
    publishDir params.out + "/filtered_indels", mode:'symlink'
    
    input:
    tuple val(pair_id), path(raw_indels)
    val(round) 

    output:
    tuple val(pair_id), \
	  path("${pair_id}_filtered_indels_${round}.vcf"), \
	  path("${pair_id}_filtered_indels_${round}.vcf.idx"), emit: filtered_indels

    script:
    """
     ${params.gatk} VariantFiltration \
        -R ${params.ref} \
        -V ${raw_indels} \
        -O ${pair_id}_filtered_indels_${round}.vcf \
	-filter-name "QD_filter" -filter "QD < 2.0" \
	-filter-name "FS_filter" -filter "FS > 200.0" \
	-filter-name "SOR_filter" -filter "SOR > 10.0"
    """
}

process splitNCigarReads{
    cache 'lenient'
    publishDir params.out + "/filtered_indels", mode:'symlink'

    input:
    tuple val(pair_id), path(raw_indels)

    output:
    tuple val(pair_id), path("${pair_id}_output.bam"), emit: split_bam
    //path("./*")
    
    script:
    """
    ${params.gatk} SplitNCigarReads \
      -R ${params.ref} \
      -I ${raw_indels} \
      -O ${pair_id}_output.bam 
    """
}
