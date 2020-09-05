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

allowed_params = ["input_vcf_ref", "input_phasefile",  "chr", "from_bp", "to_bp", "extract", "keep", "output_dir", "maf", "output", "big_time", "fasta_file"]
infosys=["scripts"]
infosoft=["bin_vcftools", "bin_vcftools", "big_time", "bin_tabix", "bin_crossmap", "bin_bcftools", "bin_beagle", "bin_shapeit"]
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
params.otheropt_beagle = ""
/*"input_col_ref", "keep", "chr", "from_bp", "to_bp", "extract", "keep", "output_dir"*/
params.input_col_ref=""
params.input_vcf_ref=""
params.keep=""
params.chr=""
params.from_bp=""
params.to_bp=""
params.extract=""
params.maf=""
params.exclude_bed=""
params.bed=""
params.vcf=""

if(params.vcf==""){
  println "This workflow requires vcf file to imputed : --vcf"
  exit 1
}

file_vcf_tofilter=Channel.fromPath(params.vcf)
vcfrefisgz=file_vcf_tofilter.getExtension()=='gz'

if(params.keep!="" || params.chr!="" || (params.chr!="" && params.to_bp!="" && params.from_bp!="") || params.exclude_bed!="" || params.bed!="" || params.){


  if(params.keep!="" || params.chr!="" || (params.chr!="" && params.to_bp!="" && params.from_bp!="") || params.exclude_bed!="" || params.bed!=""){
     if(params.keep!=""){
      keep_file_ch=Channel.fromPath(params.keep)
     }else{
      keep_file_ch=file('No_ind')
     }
     if(params.exclude_bed!=""){
       exclude_bed_ch=Channel.fromPath(params.exclude_bed)
     }else{
       exclude_bed_ch=file("exclude_bed_no")
     }
     if (params.bed=="") bed_ch=Channel.fromPath(params.bed)
     else bed_ch=file("bed_no")

     process FiltersVcfI{
      memory params.memory_vcftools
      time params.big_time
      input:
         file(file_vcf) from file_vcf_tofilter
         file(file_ind) from keep_file_ch
         file(exclude_bed) from exclude_bed_ch
         file(bedf) from bed_ch
      publishDir "${params.output_dir}/filters_imp/", overwrite:true, mode:'copy'
      output :
         file("${out_file}.recode.vcf.gz") into file_vcf_filter_1
      script :
         gvcf= vcfrefisgz  ? "--gzvcf" : "--vcf"
         chro=params.chr!=""? "--chr ${params.chr}" : ""
         end=params.to_bp!=""? "--to-bp  ${params.to_bp}" : ""
         begin=params.from_bp!=""? "--from-bp  ${params.from_bp}" : ""
         maf=params.maf!="" ? "--maf ${params.maf} " : ""
         keep=params.keep!="" ? " --keep ${file_ind} " : ""
         excl_bed=params.exclude_bed!="" ? " --exclude_bed $exclude_bed " : ""
         bed=params.bed!="" ? " --bed $bedf " : ""
         out_file="${params.output}_filt1"
         """
         ${params.bin_vcftools} $gvcf $file_vcf $chro $end $begin $maf $keep --out ${out_file} --recode --recode-INFO-all $excl_bed $bed
         ${params.bcftools} sort ${out_file}".recode.vcf"  -Oz -o  ${out_file}".recode.vcf.gz"
         """
     }

  }
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
    file_vcf_filter = file_vcf_filter_1
  }

}else{
if(vcfrefisgz){
  process gzipvcf{
      memory params.memory_vcftools
      time params.big_time
      input :
         file(file_vcf) from file_vcf_tofilter
      output :
         file("${out_file}.recode.vcf.gz") into file_vcf_filter
      script :
         out_file="${params.output}_filt1"
         """
         ${params.bcftools} sort ${file_vcf}".recode.vcf"  -Oz -o  ${out_file}".recode.vcf.gz" 
         """
  }
}else{
  file_vcf_filter=file_vcf_tofilter
}
}


