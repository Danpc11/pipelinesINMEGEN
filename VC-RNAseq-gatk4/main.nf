#!/usr/bin/env nextflow
// Flujo de trabajo llamado de variantes (GATK haplotypeCaller)
// INMEGEN 2022 
// Recuerda editar el archivo nextflow.config con los datos de tu análisis

nextflow.enable.dsl=2

include { star as align;
          markDuplicatesSpark;
          getMetrics;
          splitNCigarReads;
          bqsr;
          analyzeCovariates;
          annovarS as annovar } from "../modules.nf"
          
include { qualitycontrol;
          bootstrapping;
          bootstrapping as variantcalling } from "./subworkflow/subworkflow.nf"          

// Imprimir la ruta de algunos directorios importantes
println " "
println "Flujo de trabajo indentificación de variantes de RNA-seq"
println " "
println "Datos crudos: $params.reads"
println "Ref index: $params.ref"
println "Directorio de salida: $params.out"
println " "

workflow {
   // subflujo de trabajo para el análisis de calidad de las muestras
   qualitycontrol()
  
  // Sección para el preprocesamiento de los datos 
   Channel.fromFilePairs( "${params.reads}" )
          .tap { read_pairs_ch }
   
   align(read_pairs_ch)

   markDuplicatesSpark(align.out.aligned_reads_ch)
 
   getMetrics(markDuplicatesSpark.out.bam_for_variant_calling) 

   splitNCigarReads(markDuplicatesSpark.out.bam_for_variant_calling) 
   
   //  subflujo de trabajo para aplicar BQSR    
   round = 1
   bootstrapping(splitNCigarReads.out.split_bam,round)
   
   filteredvcf_out = bootstrapping.out.filterSnps_out.join(bootstrapping.out.filterIndels_out)
   bqsr_inlet= markDuplicatesSpark.out.bam_for_variant_calling.join(filteredvcf_out)
   bqsr_inlet.view()

   bqsr(bqsr_inlet)
   analyzeCovariates(bqsr.out.analyze_covariates)
   
   // subflujo de trabajo para el llamado de variantes
   round_f = round + 1
   variantcalling(bqsr.out.recalibrated_bam,round_f)
      
   // proceso para la anotación de variantes
   annovar(variantcalling.out.filterSnps_out)  
}
