#!/usr/bin/env nextflow
// Pipeline PureCN

nextflow.enable.dsl=2

include { coverageNormal;
          coverageTumor;
          normalDB;
          pureCN } from "../modules.nf"

// Imprimir la ruta de algunos directorios importantes
println " "
println "Directorios: PureCN pipeline"
println " "
println "Datos de archivos bam normales: $params.normal_bam"
println "Datos de archivos bam tumor: $params.tumor_bam"
println "Datos de vcfs de tumor: $params.tumor_vcfs"
println "Directorio de salida: $params.out"
println " "

workflow {

     // calculo de la covertura de muestras normales
     Channel.fromPath("${params.normal_bam}" )
          .splitCsv(sep:",", header: true)
          .map { row ->  def sample_ID = "${row.Sample}"
                         def bam_t = file("${row.Path}")
                 return [ sample_ID, bam_t]
         }.set { normal_bam_ch }

   coverageNormal(normal_bam_ch)
     
      coverage_normal=coverageNormal.out.conv_loss.toList()
   
   // Base de datos de normales
   normalDB(coverage_normal)
   
   // calculo de la covertura de las muestras de tumor
   Channel.fromPath("${params.tumor_bam}" )
          .splitCsv(sep:",", header: true)
          .map { row ->  def sample_ID = "${row.Sample}"
                         def bam_t = file("${row.Path}")
                 return [ sample_ID, bam_t]
         }.set { tumor_bam_ch }
   
   coverageTumor(tumor_bam_ch)
   
   coverage_tumor=coverageTumor.out.conv_loss.collect().flatten().collate( 2 )
   
   // calculo de CVN con pureCN
      Channel.fromPath("${params.tumor_vcfs}" )
          .splitCsv(sep:",", header: true)
          .map { row ->  def sample_ID = "${row.Sample}"
                         def vcf = "${row.Path}"
                         def stats = "${row.Stats}"
                 return [ sample_ID, vcf, stats ]
         }.set { tumor_vcf_ch}
         
    tumor_vcf_ch.join(coverage_tumor).groupTuple().flatten().collate( 4 ).set{tumor_files}
    //tumor_files.view()
         
    pureCN(tumor_files,normalDB.out.normal_files)
}
