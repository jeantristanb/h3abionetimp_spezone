#!/usr/bin/env nextflow
/*
 * Authors       :
 *
 *
 *      Jean-Tristan Brandenburg
 *
 *  2015-2020
 *
 *
 * Description  : Nextflow pipeline for imputation test.
 *
 */

//---- General definitions --------------------------------------------------//

import java.nio.file.Paths



def helps = [ 'help' : 'help' ]

allowed_params = ["input_vcf_ref", "input_col_ref", "keep", "chr", "from_bp", "to_bp", "extract", "keep", "output_dir", "maf"]

params.each { parm ->
  if (! allowed_params.contains(parm.key)) {
    println "\nUnknown parameter : Check parameter <$parm>\n";
  }
}


def params_help = new LinkedHashMap(helps)


params.queue      = 'batch'
params.work_dir   = "$HOME/"
params.output_dir = "${params.work_dir}/output"
/*"input_col_ref", "keep", "chr", "from_bp", "to_bp", "extract", "keep", "output_dir"*/
params.input_col_ref=""
params.keep=""
params.chr=""
params.from_bp=""
params.to_bp=""
params.extract=""
params.maf=""
params.bin_vcftools="vcftools"
params.memory_vcftools="10GB"

params.input_vcf_ref=""

if (params.help) {
    params.each {
    entry ->
      print "Parameter: <$entry.key>    \t Default: $entry.value"
      if (entry.key == 'help')
          println ""
      else {
        help = params_help.get(entry.key)
        if (help)
          print "\n    $help"
        println ""
      }
  }
  System.exit(-1)
}

if(params.input_vcf_ref==""){
  println "This workflow requires input_plk_toimp or input_vcf_toimp"
  exit 1
}

/*Step : selection and merge of vcf*/


//case where need to do some filtering or added information
if(params.input_col_ref!="" || params.keep!="" || params.chr!="" || (params.chr!="" && params.to_bp!="" && params.from_bp!="") || params.keep!="" || params.maf!=""){

  if(params.keep!="" || params.chr!="" || (params.chr!="" && params.to_bp!="" && params.from_bp!="")){
    file_vcf_tofilter=Channel.fromPath(params.input_vcf_ref)
    if(params.keep!=""){
    keep_file_ch=Channel.fromPath(params.keep)
    }else{
    keep_file_ch=file('No_ind')
    }
    process FiltersVcfI{
      memory params.memory_vcftools
      input: 
         file(file_vcf) from file_vcf_tofilter
         file(file_ind) from keep_file_ch 
      output :
         file(out_file) into file_vcf_filter_1
      script :
         gvcf=file_vcf.getExtension()==".gz"? "--gzvcf" : "--vcf" 
         chro=params.chr!=""? "--chr ${params.chr}" : ""
         end=params.to_bp!=""? "--to-bp  ${params.to_bp}" : "" 
         begin=params.from_bp!=""? "--from-bp  ${params.from_bp}" : "" 
         maf=params.maf!="" ? "--maf ${params.maf} " : ""
         keep=params.keep!="" ? " --keep ${file_ind} " : ""
         out_file="${params.out}_filt1"
         """
         vctools $gvf $file_vcf $chro $end $begin $maf $keep --out ${out_file}
         """ 
    } 
  } 


}else{
  file_vcf_ref_filt=Channel.fromPath(params.input_vcf_ref)
}


