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

        It was written in C++ and some libraries used might be difficut to compile. In this pipeline, I used the pre-complied binary code. 

     2. Please use the fullpath when running easyChain
        /full/path/to/easyChain/src/easyChain

     If you have any problems, please contact
 
         Zemin Ning ( zn1@sanger.ac.uk )  
