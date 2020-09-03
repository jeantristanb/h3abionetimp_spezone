# Pipeline to test different imputation with different different type of dataspecific  
## prepare ref panel : `build_ref.nf`
build reference file for XXX, XX for imputation
example command line : 
nextflow build_ref.nf --input_vcf_ref file_vcf

### arguments 
* `input_vcf_ref` :  file vcf to build panel
* filter file :
 * `input_col_ref` : file contains position to add to reference [Optional]
 * `keep` : individual to keep  [Optional]
 * `chr` : chromosome specific to keep [Optional]
 * `to_bp` : position begin [Optional]
 * `from_bp` : position begin [Optional]
 * `maf` : minor allele frequencie [Optional]
* to do :
 * `extract` : not implemented
