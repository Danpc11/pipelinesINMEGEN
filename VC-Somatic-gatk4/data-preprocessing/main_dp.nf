#!/usr/bin/env nextflow
// Pipeline preprocesamiento de datos para llamado de variantes 

nextflow.enable.dsl=2

include { fastq_to_sam;
          mark_duplicates;
          sam_to_fastq;
          align;
          merge_bam_alignment;
          mark_duplicates_spark;
          recal_bqsr } from "../../modules.nf"        

// Imprimir la ruta de algunos directorios importantes
println " "
println "Directorios: Preprocesamiento de datos para llamado de variantes"
println " "
println "Datos crudos: $params.reads"
println "Ref index: $params.ref"
println "Directorio de salida: $params.out"
println " "

workflow {
 
     //read_pairs_ch = Channel.fromFilePairs( "${params.reads}" )
    
   Channel.fromPath("${params.reads}" )
          .splitCsv(sep:"	", header: true)
          .map { row ->  def sampleID = "${row.Sample}"
                         def RG = "${row.RG}"
                         def PU = "${row.PU}"
                         def read1 = file("${row.R1}")
                         def read2 = file("${row.R2}")
                 return [ sampleID, RG, PU, read1, read2 ]
         }.set { read_pairs_ch}

  //read_pairs_ch.view()
   
   fastq_to_sam(read_pairs_ch)
       // canal que junta la salida de fastq to sam
       ch_fqtsam=fastq_to_sam.out.fastq_to_sam_ch.collect().flatten().collate( 2 )
 
   mark_duplicates(fastq_to_sam.out.fastq_to_sam_ch)
   
   sam_to_fastq(mark_duplicates.out.mark_duplicates_bam)
          
   align(sam_to_fastq.out.sam_to_fastq_ch)
       // canales que juntan las salidas de fastq to sam y align para hacer merge por nombre de muestra
       // los acomoda en un formato [sampleID, path fastqtosam, path align bam]
       ch_aln=align.out.aligned_reads.collect().flatten().collate( 2 )
       ch_aln.join(ch_fqtsam).groupTuple().flatten().collate( 3 ).set{joinfmerge}
   
      merge_bam_alignment(joinfmerge)
   
   mark_duplicates_spark(merge_bam_alignment.out.bam_alignment_ch)

   recal_bqsr(mark_duplicates_spark.out.bam_for_variant_calling)
   
}
