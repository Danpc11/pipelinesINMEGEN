# Flujo de trabajo llamado de variantes de línea germinal utilizando NextFlow y GATK4 (pipeline VC_Germline)

Este pipeline realiza la indentificación conjunta de variantes a partir de archivos de secuenciación masiva (WGS/WES) en un formato *fastq*.

#### Las herramientas necesarias para correr este flujo de trabajo son:
> 
> - NextFlow (22.04.5)
> - R (4.2.0) 
> - FastQC (0.11.9) 
> - MultiQC (1.11)
> - Openjdk (11.0.13 o superior)
> - GATK (4.2.6.1)
> - BWA (0.7.17-r 1188)
> - Picard Tools (2.0.1)
> - Samtools (1.6.0)
> - SnpEff (5.1)
> - Annovar
> - bcftools (1.14)
> 

Opcional, contar con un archivo csv o txt con la información y el ReadGroup de las muestras. 
Los archivos como el índice de [BWA](http://bio-bwa.sourceforge.net/) y los archivos de recalibración de BQSR y VQSR se pueden descargar del [bundle de GATK](https://console.cloud.google.com/storage/browser/genomics-public-data/resources/broad/hg38/v0;tab=objects?prefix=&forceOnObjectsSortingFiltering=false).  

## Instructivo de uso

### Para Ejecutar por primera vez el pipeline 
 
 1. Verificar que el entorno de trabajo cuente con las herramientas descritas en la sección herramientas necesarias, en caso de que no se tengan instaladas, se recomienda generar un microambiente conda de NextFlow e instalar ahí las herramientas necesarias (excepto, GATK, Picard y Snpeff, de estas se especifica la ruta en el archivo de configuración). Para realizar lo anterior se recomienda seguir los siguientes pasos:	
 
			conda install -c bioconda nextflow 
			conda install -c bioconda bwa
			conda install -c bioconda samtools
			conda update samtools
			conda install -c bioconda bcftools 

 2. Descargar GATK (4.2.6.1), Picard (2.0.1), Annovar y SnpEff (5.1)
 3. Verificar que los archivos de nextflow estén disponibles (descargar de este repositorio), los archivos necesarios son:
	- archivo: main.nf
	- archivo: modules.nf
	- archivo: subworkflow.nf
	- archivo: nextflow.config
 4. Verificar la ruta del directorio de Openjdk 8 (ver sección microambiente conda para GATK)
 5. Archivos a analizar en un formato *fastq*
 6. La ruta con los índices de BWA y SamTools
 7. La ruta con los archivos de recalibración de BQSR y VQSR
 8. La ruta con las bases de datos de annovar

Una vez que se cuenta con lo enunciado en las secciones anteriores se puede proceder a correr el pipeline VC_GATK siguiendo los siguientes pasos:

### Ejecutar flujo de trabajo VC_GATK

 1. Generar un directorio en el que se concentrará todo lo necesario para el análisis, se recomienda un nombre que incluya la fecha del análisis y una breve descripción del diseño experimental, p.ej. AAAA-M-D_DEA-nombre del ensayo. Esta carpeta se denominará directorio principal
 2. Opcional: Si los datos curdos (archivos *fastq*) se encuentran en diferentes directorios, moverlos a una sóla carpeta utilizando
 
			find <carpeta_origen> -name "*fastq.gz" -print0 | xargs -0 -I {} mv {} <carpeta_destino>
			
 3. Editar el archivo de nexflow.config con la siguiente información:
	- Ruta de los archivos *fastq*
	- Ruta del directorio de salida de nextflow
	- Nombre del proyecto 
	- Ruta del índice de BWA (*.idx)
	- Ruta del ejecutable de GATK
	- Ruta del ejecutable de Picard
	- Ruta del ejecutable de SnoEff
	- Ruta de los vcfs de BQSR y VQSR
	- Ruta de las bases de datos de annovar
 4. En el directorio donde están alojados los archivos de nexflow ejecutar el comando:
 
			nextflow run <pipeline name>
			
Si se quiere un informe del rendimiento de Nextflow se puede optar por:
 
			nextflow run <pipeline name> -with-report <file name>
			
Una vez ejecutado el pipeline, verificar que la carpeta de resultados dentro del directorio de salida contenga los archivos *.vcf y los informes *.html.
 
**NOTA 1:** En caso de ocurrir un error debido a una ruta mal escrita o algún evento similar, se puede editar los archivos de nextflow corrigiendo los errores, guardarlo de nuevo y resumir el pipeline del análisis utilizando la opción -resume, como se muestra a continuación: 

			nextflow run <pipeline name> -resume

### Microambiente conda para GATK 

En este microambiente se pueden correr las diferentes herramientas de GATK que utilizan la versión 8 de Openjdk. Además, en este microambiente se puede correr el pipeline de [Mohammed Khalfan](https://github.com/gencorefacility/variant-calling-pipeline-gatk4). Sin embargo, es un ambiente alterno al que se utiliza para correr el pipeline VC_GATK

Los pasos para montar este microambiente son:
1. Instalar GATK siguiendo las instrucciones de este [link](https://gatk.broadinstitute.org/hc/en-us/articles/360035889851--How-to-Install-and-use-Conda-for-GATK4)
2. Activar el ambiente conda GATK con el siguiente comando

			conda activate GATK
			
3. Instalar la correcta versión de Java 8, utilizando el siguiente comando
			
			conda install openjdk==8.0.332=h166bdaf_0 

 *NOTA:* En caso de error intentar realizar la siguiente configuración de conda y volver a intentar la instalación 

			conda config --add channels defaults
			conda config --add channels bioconda
			conda config --add channels conda-forge
			
4. Utilizar conda para instalar BWA, Samtools y Nextflow (siguiendo ese orden) con los comandos

			conda install -c bioconda bwa
			conda install -c bioconda samtools
			conda install -c bioconda nextflow

5. Utilizar la versión v4_3t de SnpEff

Como se puede ver en el orden de los pasos, en este microambiente NextFlow es la última herramienta que se instala (por cuestiones de compatibilidad), mientras que en el microambiente del pipeline es la primera
