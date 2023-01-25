include { fastqc;
          multiqc;
          haplotypeCaller_vcr as haplotypeCaller;
          selectVariants_vcr as selectVariants;
          filterSnps_vcr as filterSnps;
          filterIndels_vcr as filterIndels} from "../../modules.nf"

// Sub flujo para el analisis de calidad de las muestras
workflow qualitycontrol {

   data_fq = Channel.fromFilePairs("${params.reads}")
                    .ifEmpty { error "Cannot find any reads matching: ${params.reads}"  }
                       
   fastqc(data_fq)
    
   multiqc(fastqc.out.fq_files.collect())   
}

// Sub flujo de trabajo para recalibración y el descubrimiento de variantes
workflow bootstrapping {

   take: data_1 
         round_1
   
   main:
           
   haplotypeCaller(data_1,round_1)
   
   selectVariants(haplotypeCaller.out.hc_output,round_1)
   
   filterSnps(selectVariants.out.raw_snps_ch,round_1)
   
   filterIndels(selectVariants.out.raw_indels_ch,round_1)
         
   emit:
   round_1 = round_1
   filterSnps_out = filterSnps.out.filtered_snps
   filterIndels_out = filterIndels.out.filtered_indels  
   //filteredvcf_out = filterSnps.out.filtered_snps.join(filterIndels.out.filtered_indels) 
}
