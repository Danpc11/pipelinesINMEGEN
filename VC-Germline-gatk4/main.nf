#!/usr/bin/env nextflow
// Flujo de trabajo indentificación conjunta de variantes con GATK4
// INMEGEN 2022
// Versión 1 (sin bootstraping)

nextflow.enable.dsl=2

include { alignG;
          mergeSam;
          markDuplicatesSpark;
          getMetrics;
          recal_bqsr;
          analyzeCovariates;
          haplotypeCallerERC;
          genomicsDBimport;
          genotypeGVCFs;
          selectVariantsConjunto;
          VQSRsnps;
          VQSRindels;
          snpEff;
          annovarG;
          splitVCFs as splitVCF_snpeff;
          splitVCFs as splitVCF_annovar  } from "../modules.nf"

include { qualitycontrol } from "./subworkflow/subworkflow.nf"

// Imprimir la ruta de algunos directorios importantes
println " "
println "Directorios del flujo de trabajo: Indetificación conjunta de variantes con GATK"
println "Tipo: Germinal"
println " "
println "Datos crudos: $params.reads"
println "Ref index: $params.ref"
println "Directorio de salida: $params.out"
println "Nombre del proyecto: $params.project_name"
println " "

workflow {
   // subflujo de trabajo para el análisis de calidad de las muestras
   qualitycontrol()

  // Sección para el preprocesamiento de los datos
     num_samples = 0
   Channel.fromFilePairs( "${params.reads}" )
          .tap { read_pairs_ch }
          .subscribe({ num_samples += 1 })

   alignG(read_pairs_ch)

     xa = alignG.out.aligned_reads_ch.collect().flatten().collate( 2 )
     xa.map { a , b -> def key = a.toString().tokenize('_').get(1)
                        return tuple(key, b)
            }.groupTuple() | mergeSam

   markDuplicatesSpark(mergeSam.out.merged_sam_ch)

   getMetrics(markDuplicatesSpark.out.bam_for_variant_calling)

   recal_bqsr(markDuplicatesSpark.out.bam_for_variant_calling)

   analyzeCovariates(recal_bqsr.out.analyze_covariates)

   // Indentificación conjunta de variates
   haplotypeCallerERC(recal_bqsr.out.recalibrated_bam)

     hc_files = haplotypeCallerERC.out.hc_erc_out.toList()
     project_id="${params.project_name}"

   genomicsDBimport(hc_files,project_id)

   genotypeGVCFs(genomicsDBimport.out.genomics_db)

   selectVariantsConjunto(genotypeGVCFs.out.gvcfs_out)

   VQSRsnps(selectVariantsConjunto.out.snps_ch)

   VQSRindels(selectVariantsConjunto.out.indels_ch)

   // Anotación de variantes
   snpEff(VQSRsnps.out.snps_filt_ch,VQSRindels.out.indels_filt_ch)
   splitVCF_snpeff(snpEff.out.snpeff_ch)

   annovarG(VQSRsnps.out.snps_filt_ch,VQSRindels.out.indels_filt_ch)
   splitVCF_annovar(annovarG.out.annovar_ch)
} 
