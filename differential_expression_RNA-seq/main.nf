#!/usr/bin/env nextflow

// Flujo de trabajo Cuantificación y análisis de expresión diferencial
// INMEGEN 2022 

nextflow.enable.dsl=2

include {fastqc;
	 trim_Galore;
	 multiqc_DEA as multiqc;
         kallisto;
         tximport_deseq2;
         tximport_q } from "../modules.nf"
         
// Imprimir algunos directorios importantes
println " "
println "Cuantificación y Análisis de Expresión"
println " "
println "Ruta de los datos crudos: ${params.datadir}"
println "Realizar análsis de expresión diferencial: ${params.QDEA}"
println "Ruta del directorio de salida: ${params.out}"
println " "

workflow {

    data_fq = Channel.fromFilePairs("${params.datadir}")
                     .ifEmpty { error "Cannot find any reads matching: ${params.reads}"  } 
       
    fastqc(data_fq)
    
    trim_Galore(data_fq)
       
    multiqc(fastqc.out.fq_files.collect(),trim_Galore.out.trim_html.collect())
    
    kallisto(trim_Galore.out.trim_fq)
    
    if ("${params.QDEA}" == true){ 
    tximport_deseq2(kallisto.out.abundance_h5.collect()) 
    } 
    else { 
    tximport_q(kallisto.out.abundance_h5.collect()) 
    }
}
