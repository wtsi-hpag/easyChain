# easyChain 
Pipeline to generate chain file for assembly coordinate conversion.

To perform whole genome alignment between target assembly_1.fasta (GRCh37-assembly) and reference assembly_2.fasta (GRCh38-assembly) following the steps:

           1. The target assembly is shredded into chunks of 20000 bases 
           2. The 20000 bases chunks are mapped against the reference assembly

### Download and Compile:
Requirements for compiling: gcc

		$ git clone https://github.com/SangerHpag/easyChain.git
		$ cd easyChain 
		$ ./install.sh
		
If everything compiled saccessfully you must see the final comment: 
		"Congrats: installation successful!"		

(Tested with gcc-4.9.2, gcc-4.9.4, gcc-4.8.1, gcc-6.0.2) 

#### External packages
The genome aligner BWA (http://bio-bwa.sourceforge.net) and SMALT (http://www.sanger.ac.uk/science/tools/smalt-0) are downloaded and compiled by easyChain.

#### Run:

           $ /full/path/to/easyChain/src/easyChain -nodes <nodes> -shred <shred_length> \
	   	      </full/path/to/assembly_1.fasta> </full/path/to/assembly_2.fasta> \ 
		      <shred-align.out>
           
           where:
	          /full/path/to/assembly_1.fasta: full path to the assembly file to be considered as "GRCh37 assembly"
	     	  /full/path/to/assembly_2.fasta:  full path to the assembly file to be considered as "GRCh38 assembly"
	     	  shred-align.out:   output name for the genome wide alignment. 
	     
	       parameters:
             nodes:    number of CPUs requested  [ default = 30 ]
             shred:    length of shredded fragments [ default = 20000 ]
             
 
