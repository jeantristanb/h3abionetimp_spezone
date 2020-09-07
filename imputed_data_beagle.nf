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
params.output_dir = "${params.work_dir}/imputed_file"
params.output = "imputed"
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
params.prephase=1

params.bin_beagle="beagle"
params.bin_tabix="tabix"
params.bin_vcftools="vcftools"
params.bin_bcftools="bcftools"
params.bin_shapeit="shapeit"
params.bin_eagle="eagle"

params.memory_vcftools="10GB"
params.memory_tabix="10GB"
params.cpus_other=4
params.bin_crossmap="~/.local/bin/CrossMap_beta.py"

params.genetic_map=""
params.genetic_map_beagle=""

params.vcf_ref_norm=""


if(params.vcf==""){
  println "This workflow requires vcf file to imputed : --vcf"
  exit 1
}

println params.vcf
file_vcf_tofilter=Channel.fromPath(params.vcf)
vcfrefisgz=file(params.vcf).getExtension()=='gz'
//vcfrefisgz=params.vcf.getExtension()=='gz'


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
     if (params.bed!="") bed_ch=Channel.fromPath(params.bed)
     else bed_ch=file("bed_no")

     process FiltersVcfI{
      memory params.memory_vcftools
      time params.big_time
      input:
         file(file_vcf) from file_vcf_tofilter
         file(file_ind) from keep_file_ch
         file(exclude_bed) from exclude_bed_ch
         file(bedf) from bed_ch
      publishDir "${params.output_dir}/vcffilter/", overwrite:true, mode:'copy'
      output :
         set file("${out_file}.recode.vcf.gz"),  file("${out_file}.recode.vcf.gz.csi") into file_vcf_filter_1
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
         ${params.bin_bcftools} sort ${out_file}".recode.vcf"  -Oz -o  ${out_file}".recode.vcf.gz"
         ${params.bin_bcftools} index ${out_file}".recode.vcf.gz"
         """
     }

}else{
 if(vcfrefisgz){
  process gzipvcf{
      memory params.memory_vcftools
      time params.big_time
      input :
         file(file_vcf) from file_vcf_tofilter
      output :
         set file("${out_file}.recode.vcf.gz"),file("${out_file}.recode.vcf.gz.csi") into file_vcf_filter
      script :
         out_file="${params.output}_filt1"
         """
         ${params.bin_bcftools} sort ${file_vcf}".recode.vcf"  -Oz -o  ${out_file}".recode.vcf.gz" 
         """
  }
}else{
  file_vcf_filter=file_vcf_tofilter
}
}

/*phaseing inital data*/
if(params.prephase==1){


if(params.prephase_eagle==0){

else{
eagle_genetic_map_ch=Channel.fromPath(params.genetic_map)
ref_vcf_ch = Channel.fromPath(params.vcf_ref_norm)
process prephase_eagle{
    input : 
      file(file_vcf) from file_vcf_filter
      file(eagle_genetic_map) from eagle_genetic_map_ch
      file(ref_vcf) from ref_vcf_ch
    output :
      file(fileout) into file_vcf_prephase
    script :
        chro=params.chr!=""? "--chrom=${params.chr}" : ""
        begin=params.to_bp!=""? "--bpStart=${params.to_bp}" : ""
        end=params.from_bp!=""? "--bpEnd=${params.from_bp}" : ""
        file_out=${params.output"_prephase_eagle"}
        """
        ${params.bin_eagle}\
                --vcfTarget=${file_vcf} \
                --geneticMapFile=${eagle_genetic_map} \
                --vcfRef=${ref_vcf} \
                --vcfOutFormat=z \
                --noImpMissing \
                $chro \
                $begin \
                $end \ 
                --bpFlanking=${params.buffer_size} \
                --outPrefix=${file_out} 2>&1 | tee ${file_out}.log
       """
}
}
}else{
file_vcf_prephase=file_vcf_filter
}
/*imputed */
if(params.refpanel_bref3!=""){
}



