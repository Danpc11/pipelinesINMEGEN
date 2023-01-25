#!/usr/bin/env nextflow
// Pipeline indentificación de variantes somaticas (Mutec2)

nextflow.enable.dsl=2

include { mutect2;
          calculateContamination;
          filterMutectCalls } from "../modules.nf"

// Imprimir la ruta de algunos directorios importantes
println " "
println "Directorios: Indentificación de variantes(somaticas)"
println " "
println "Datos crudos: $params.reads_tumor"
println "Ref index: $params.ref"
println "Panel de normales: $params.panelnormales"
println "Directorio de salida: $params.out"
println " "

workflow {

     Channel.fromPath("${params.reads_tumor}" )
          .splitCsv(sep:",", header: true)
          .map { row ->  def sample_ID = "${row.Sample}"
                         def bam_t = file("${row.Path}")
                 return [ sample_ID, bam_t]
         }.set { ready_tumor_bam_ch }

   mutect2(ready_tumor_bam_ch)

     unfilt=mutect2.out.unfilt_vcf.collect().flatten().collate( 3 )
     //unfilt.view()

   calculateContamination(ready_tumor_bam_ch)

      contables=calculateContamination.out.cont_tables.collect().flatten().collate( 3 )
      //contables.view()

   unfilt.join(contables).groupTuple().flatten().collate( 5 ).set{forfilter}
    //forfilter.view()

   filterMutectCalls(forfilter)
}
