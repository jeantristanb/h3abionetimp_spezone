# Pipeline to test different imputation with different different type of dataspecific  
## prepare ref panel : `build_ref.nf`
build reference file for XXX, XX for imputation
example command line : 
nextflow build_ref.nf --input_vcf_ref file_vcf

### arguments  
* filter file :
 * `input_col_ref` : file contains position to add to reference [Optional]
  * first column must be same identifiant that vcf, each column must be rs\_chr\_pos\_1 for allele 1 and rs\_chr\_pos\_2 for allele 2
  * used R to merge both
 * `keep` : individual to keep  [Optional]
 * `chr` : chromosome specific to keep [Optional]
 * `to_bp` : position begin [Optional]
 * `from_bp` : position begin [Optional]
 * `maf` : minor allele frequencie [Optional]
 * `exclude_bed`
 * `bed` : 
* to do :
 * `extract` : not implemented
* convert between format (optional):
 * Used crossmap in python to convert between 2 references hg38 => hg19, hg19 => hg38 
 * `convert_file` : file for crossmap for instance [file uscs](http://hgdownload.cse.ucsc.edu/goldenpath/hg38/liftOver/) [default NULL] if file not present will be not convert
 * `fasta_file` : fasta file from final reference 
 * `bin_crossmap` : [default : ~/.local/bin/CrossMap_beta.py]
General Arguements for imputation of reference panel :
 * `genetic_map` : genetic map, 
* beagle :
 * `bin_beagle` : binary beagle [default beagle]
 * `otheropt_beagle` : other option 
* shapeit :
 * `thr_shapeit` : threshold for shapeit [default 0.95]
 * `bin_shapeit`

## imputed vcf : `build_ref.nf`
## need
* general : nextflow,  
* `build_ref.nf` :  tabix, vcftools, CrossMapp
