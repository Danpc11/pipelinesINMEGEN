include { fastqc;
          multiqc;
          haplotypeCaller;
          selectVariants;
          filterSnps;
          filterIndels } from "../../modules.nf"

// Sub flujo para el analisis de calidad de las muestras
workflow qualitycontrol {

   data_fq = Channel.fromFilePairs("${params.reads}")
                    .ifEmpty { error "Cannot find any reads matching: ${params.reads}"  }
                       
   fastqc(data_fq)
    
   multiqc(fastqc.out.fq_files.collect())   
}

// Sub flujo de trabajo para recalibración (bootstrapping)
workflow bootstrapping {

   take: data_1 
   
   main:
           
   haplotypeCaller(data_1)
   
   selectVariants(haplotypeCaller.out.hc_output)
   
   filterSnps(selectVariants.out.raw_snps_ch)
   
   filterIndels(selectVariants.out.raw_indels_ch)
         
   emit:
   filterSnps_out = filterSnps.out.filtered_snps
   filterIndels_out = filterIndels.out.filtered_indels     
} 
