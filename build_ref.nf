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

allowed_params = ["input_vcf_ref", "input_col_ref", "keep", "chr", "from_bp", "to_bp", "extract", "keep", "output_dir", "maf", "output", "big_time", "convert_file", "fasta_file"]
infosys=["scripts"]
infosoft=["bin_vcftools", "bin_vcftools", "big_time", "bin_tabix", "bin_crossmap"]
allowed_params+=infosoft

params.each { parm ->
  if (! allowed_params.contains(parm.key)) {
    println "\nUnknown parameter : Check parameter <$parm>\n";
  }
}


def params_help = new LinkedHashMap(helps)


params.queue      = 'batch'
params.work_dir   = "$HOME/"
params.output_dir = "${params.work_dir}/reference_panel"
params.output = "refpanel"
/*"input_col_ref", "keep", "chr", "from_bp", "to_bp", "extract", "keep", "output_dir"*/
params.input_col_ref=""
params.keep=""
params.chr=""
params.from_bp=""
params.to_bp=""
params.extract=""
params.maf=""
params.convert_file=""
params.bin_tabix="tabix"
params.bin_vcftools="vcftools"
params.memory_vcftools="10GB"
params.memory_tabix="10GB"
params.bin_crossmap="~/.local/bin/CrossMap_beta.py"

params.input_vcf_ref=""

if (params.help) {
    params.each {
    entry ->
      print "Parameter: <$entry.key>    \t Default: $entry.value"
      if (entry.key == 'help')
   bin_vcftools       println ""
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

/*define if ftp or not*/
vcfrefext=0
tmpref=params.input_vcf_ref.split("/")[0]
if(tmpref=="http:" || tmpref=="https:" || tmpref=="ftp:"){
      vcfrefext=1
}

/*gzip or not*/
//vcfrefisgz=params.input_vcf_ref =~ /.gz$/ 

if(vcfrefext){
 process GetDataFTP{
 memory params.memory_tabix
  output :
     file(newvcf) into file_vcf_tofilter  
  script :
    newvcf=params.output+"_sub.vcf"
    posinf=params.chr!=""? "${params.chr}" : "" 
    posinf=(params.to_bp!="" && params.from_bp!="") ? " $posinf:${params.from_bp}-${params.to_bp}" : posinf
    """ 
    ${params.bin_tabix} -h ${params.input_vcf_ref} $posinf > $newvcf
    """
 }


}else{
      file_vcf_tofilter=Channel.fromPath(params.input_vcf_ref)
}

vcfrefisgz=file_vcf_tofilter.getExtension()=='gz'
//case where need to do some filtering or added information
if(params.input_col_ref!="" || params.keep!="" || params.chr!="" || (params.chr!="" && params.to_bp!="" && params.from_bp!="") || params.keep!="" || params.maf!="" || params.convert_file!=""){

  if(params.keep!="" || params.chr!="" || (params.chr!="" && params.to_bp!="" && params.from_bp!="")){
    if(params.keep!=""){
      keep_file_ch=Channel.fromPath(params.keep)
    }else{
      keep_file_ch=file('No_ind')
    }

    process FiltersVcfI{
      memory params.memory_vcftools
      time params.big_time
      input: 
         file(file_vcf) from file_vcf_tofilter
         file(file_ind) from keep_file_ch 
      output :
         file("${out_file}.recode.vcf") into file_vcf_filter_1
      script :
         gvcf= vcfrefisgz  ? "--gzvcf" : "--vcf" 
         chro=params.chr!=""? "--chr ${params.chr}" : ""
         end=params.to_bp!=""? "--to-bp  ${params.to_bp}" : "" 
         begin=params.from_bp!=""? "--from-bp  ${params.from_bp}" : "" 
         maf=params.maf!="" ? "--maf ${params.maf} " : ""
         keep=params.keep!="" ? " --keep ${file_ind} " : ""
         out_file="${params.output}_filt1"
         """
         ${params.bin_vcftools} $gvcf $file_vcf $chro $end $begin $maf $keep --out ${out_file} --recode --recode-INFO-all
         """ 
    } 
    vcfrefisgz=0
  }else{
   file_vcf_filter_1 = file_vcf_tofilter
  } 
  /*transformation of position on other position */
  if (params.convert_file!=""){
    if(params.fasta_file==""){
      println "error fasta file not found" 
      exit 1
    }
    fastafile_ch=Channel.fromPath(params.fasta_file)
    file_conv_ch=Channel.fromPath(params.convert_file)
    process ConvertPosition{
      input :
         file(vcfI) from file_vcf_filter_1
         file(convert) from file_conv_ch
         file(fasta) from fastafile_ch
       publishDir "${params.output_dir}/crossmap_out/", overwrite:true, mode:'copy'
       output :
          file("${vcffinal}.unmap")
          file(vcffinal) into file_vcf_filter_2 
       script :
        vcffinal="${params.output}_filt1_newpos.vcf"
        """
        ${params.bin_crossmap}  vcf $convert $vcfI $fasta $vcffinal
        """
    }

  }else{
    file_vcf_filter_2 = file_vcf_filter_1
  }
  if(params.input_col_ref!=""){
     inputcolref=Channel.fromPath(params.input_col_ref)
     process AddPosition {
        input:
          file(filevcfi) from file_vcf_filter_2
          file(colreffile) from inputcolref  
        script :
          """ 
          """
     }   
  }else{
    file_vcf_filter_final=file_vcf_filter_2
  }

}else{
  file_vcf_filter_final=file_vcf_tofilter
}


