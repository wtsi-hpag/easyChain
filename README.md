# easyChain 
Pipeline to generate chain file for assembly coordinate conversion.

To perform whole genome alignment between target assembly_1.fasta (GRCh37-assembly) and reference assembly_2.fasta (GRCh38-assembly) following the steps:

           1. The target assembly is shredded into chunks of 20000 bases 
           2. The 20000 bases chunks are mapped against the reference assembly
           3. Generate a standard chain file shredOut.chain using the alignment file 
      
### Download and Compile:
Requirements for compiling: gcc

		$ git clone https://github.com/wtsi-hpag/easyChain.git
		$ cd easyChain 
		$ bash install.sh
		
If everything compiled saccessfully you must see the final comment: 
		"Congrats: installation successful!"		

(Tested with gcc-4.9.2, gcc-4.9.4, gcc-4.8.1, gcc-6.0.2) 

#### External packages
The genome aligner BWA (http://bio-bwa.sourceforge.net) and SMALT (http://www.sanger.ac.uk/science/tools/smalt-0) are downloaded and compiled by easyChain.

#### Run:

           $ /full/path/to/easyChain/src/easyChain -nodes <nodes> -shred <shred_length> \
	   	      </full/path/to/assembly_1.fasta> </full/path/to/assembly_2.fasta> <shredOut.chain>\ 
           
           where:
	          /full/path/to/assembly_1.fasta: full path to the assembly file to be considered as "GRCh37 assembly"
	     	  /full/path/to/assembly_2.fasta:  full path to the assembly file to be considered as "GRCh38 assembly"
	     	  shredOut.chain:   output name for the standard chain file. 
	     
	       parameters:
             nodes:    number of CPUs requested  [ default = 30 ]
             shred:    length of shredded fragments [ default = 20000 ]
             output:   output file (1) alignment only; (2) standard chain file [ default = 2 ]
             
#### Note
     1. The shred2chain part is developed by Yongji Liu in Beijing, China, see

        https://github.com/liu-yongji/shred2chain

        It was written in C++ and some libraries used might be difficut to compile. 
        In this pipeline, I used the pre-complied binary code. 

     2. Please use the fullpath when running easyChain
        /full/path/to/easyChain/src/easyChain


### checkError
We provide a pipeline to check conversion errors made by different chain files
Before using this tool, you need to have the VCF files ready:

     1. Download and unzip all the VCF files from the 1000 Genome project for GRCH37;( chr1~22，X,Y)
        http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/
        make a directory grch37_vcf and copy all the unzipped GRCh37 VCF files there

     2. Download and unzip all the VCF files from the 1000 Genome project for GRCH38.( chr1~22，X,Y)
        http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/supporting/GRCh38_positions/
        make a directory grch38_vcf and copy all the unzipped GRCh38 VCF files there

     3. Disease databases such as ClinVar, Gwas_catlog, hgmd and omim
        We have included these 4 datasets in the pipeline and they will be ready after installation! 
        No action is needed for users
     
     4. VCF files have to be the annotated files with RS numbers assigned to each called varrant
        Self generated VCFs without annotation will not work

#### USAGE 

	[Usage]: ./checkError [VCF_Files_folder] [reference_VCF_Files_folder] [chain_file] [output_folder] \ 
	[Example]: ./checkError grch37_vcf grch38_vcf hg19ToHg38.over.chain output_result \ 

   	grch37_vcf            - The folder which contains some or all the VCF files for GRCh37 \
   	grch38_vcf            - The folder which contains ALL the VCF files for GRCh38    \
   	hg19ToHg38.over.chain - The chain file selected                                   \
   	output_result         - The folder with the output results                        \

#### Note 

You need to create a directory for output results. If not exists, the code will generate a temp one \ 
Five files will be generated after processing for each input VCF file, they are:     \

 	[1]  xxx_SNP.bed:                   - The bed file extracted from xxx.vcf with tag "VT=SNP".   \
 	[2]  xxx_SNP_genegos.bed:           - The file after coordinate conversion.                    \
 	[3]  xxx_SNP_genegos.unmap:         - The content that could not be converted.                 \
 	[4]  xxx_SNP_genegos_error.dat:     - The file contains all the error sites.                   \
 	[5]  xxx_SNP_genegos_error_db.txt:  - The file contains all error sites in import databases.   \

#### Further information

     1. The checkError pipeline is developed by Yongji Liu in Beijing, China, see
        https://github.com/liu-yongji/checkbederror

     2. Zemin Ning integrated all the codes and disease databases here for better download and installation. 

     If you have any problems, please contact
 
         Zemin Ning ( zn1@sanger.ac.uk )  
