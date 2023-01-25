#!/usr/bin/env nextflow
// Pipeline preprocesamiento de datos para llamado de variantes 

nextflow.enable.dsl=2

include { mutect2forPanelofNormals;
          genomicsDBImport;
          createSomaticPanelofNormals } from "../../modules.nf" 

// Imprimir la ruta de algunos directorios importantes
println " "
println "Directorios: Generar panel de normales para Mutec2"
println " "
println "Datos crudos (muestras normales): $params.reads"
println "Ref index: $params.ref"
println "Directorio de salida: $params.out"
println " "

workflow {
 
   Channel.fromPath("${params.reads}" )
          .splitCsv(sep:",", header: true)
          .map { row ->  def sampleID = "${row.Sample}"
                         def bam_f = file("${row.Path}")
                 return [ sampleID, bam_f]
         }.set { ready_bam_ch }

   mutect2forPanelofNormals(ready_bam_ch)
   
     vcf_files = mutect2forPanelofNormals.out.mtf_PON_out.toList()   

   genomicsDBImport(vcf_files)
   
   createSomaticPanelofNormals(genomicsDBImport.out.lfor_PON_out)
}
