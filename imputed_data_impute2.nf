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

allowed_params = ["vcf_ref", "input_phasefile",  "chr", "from_bp", "to_bp", "extract", "keep", "output_dir", "maf", "output", "big_time", "fasta_file"]
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
params.haps_ref= ""
/*"input_col_ref", "keep", "chr", "from_bp", "to_bp", "extract", "keep", "output_dir"*/
params.keep=""
params.chr=""
params.from_bp=""
params.to_bp=""
params.extract=""
params.maf=""
params.exclude_bed=""
params.bed=""
params.vcf=""
params.vcf_ref=""
params.prephase=0
params.thr_shapeit=0.95

params.effective_size=20000
params.bin_beagle="beagle"
params.bin_tabix="tabix"
params.bin_vcftools="vcftools"
params.bin_bcftools="bcftools"
params.bin_shapeit="shapeit"
params.bin_eagle="eagle"
params.bin_bref3="bref3"
params.bin_impute2="impute2"
params.buffer_kb=100
params.chro=""

params.memory_vcftools="10GB"
params.memory_tabix="10GB"
params.cpus_other=4
params.bin_crossmap="~/.local/bin/CrossMap_beta.py"

params.genetic_map="tmp"
params.genetic_map_beagle=""

params.vcf_ref_norm=""


if(params.vcf==""){
  println "This workflow requires vcf file to imputed : --vcf"
  exit 1
}

println params.vcf
file_vcf_tofilter=Channel.fromPath(params.vcf, checkIfExists:true)
vcfrefisgz=file(params.vcf).getExtension()=='gz'
//vcfrefisgz=params.vcf.getExtension()=='gz'


if(params.keep!="" || params.chr!="" || (params.chr!="" && params.to_bp!="" && params.from_bp!="") || params.exclude_bed!="" || params.bed!=""){
     if(params.keep!=""){
      keep_file_ch=Channel.fromPath(params.keep,checkIfExists:true)
     }else{
      keep_file_ch=file('No_ind')
     }
     if(params.exclude_bed!=""){
       exclude_bed_ch=Channel.fromPath(params.exclude_bed, checkIfExists:true)
     }else{
       exclude_bed_ch=file("exclude_bed_no")
     }
     if (params.bed!="") bed_ch=Channel.fromPath(params.bed,checkIfExists:true)
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
         file("${out_file}.recode.vcf.gz") into file_vcf_filter
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
         """
     }

}else{
 if(vcfrefisgz){
  process gzipvcf{
      memory params.memory_vcftools
      time params.big_time
      input :
         set file(file_vcf) from file_vcf_tofilter 
      output :
         set file("${out_file}.recode.vcf.gz") into file_vcf_filter
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

process normgvcf{
     memory params.memory_vcftools
     time params.big_time
     input :
       file(file_vcf) from file_vcf_filter
      output :
       set file(out_file), file("${out_file}.csi") into file_vcf_norm
      script : 
       out_file="${params.output}_norm.vcf.gz"
       """
       ${params.bin_bcftools} norm -m -any $file_vcf -Oz -o $out_file
       ${params.bin_bcftools}  index $out_file
       """ 
}
/*phaseing inital data*/

genetic_map_ch=Channel.fromPath(params.genetic_map, checkIfExists:true)
 process prephase_shapeit{
  cpus params.cpus_other
  input :
      set file(file_vcf), file(vcffilecsi) from file_vcf_norm
      file(genetic_map) from genetic_map_ch
  publishDir "${params.output_dir}/shapeitprephase/", overwrite:true, mode:'copy'
  output :
      file("${fileout}.log") 
      set file("${fileout}.haps"), file("${fileout}.sample")   into file_haps_prephase
  script :
      chro=params.chr!=""? "--chrom=${params.chr}" : ""
      begin=params.to_bp!=""? "--bpStart=${params.from_bp}" : ""
      end=params.from_bp!=""? "--bpEnd=${params.to_bp}" : ""
      fileout="${params.output}_prephase_shapeit"
      """
      ${params.bin_shapeit} \
       -V $file_vcf\
       --input-map $genetic_map \
       --input-thr ${params.thr_shapeit} \
       --output-max $fileout".haps" $fileout".sample" \
       --thread ${params.cpus_other} \
       --effective-size ${params.effective_size} \
       --output-log ${fileout}".log" \
       --force
      """ 
 }

/*imputed */
genetic_map_ch_2=Channel.fromPath(params.genetic_map, checkIfExists:true)
if(params.haps_ref!=""){

haps_ref_file=Channel.fromPath("${params.haps_ref}.haps", checkIfExists:true)
leg_ref_file=Channel.fromPath("${params.haps_ref}.leg", checkIfExists:true)
}else{
vcf_ref=Channel.fromPath(params.vcf_ref, checkIfExists:true)
 process convertref{
  input :
   file(vcf) from vcf_ref
  output :
   file("${headout}.hap.gz") into haps_ref_file
   file("${headout}.legend.gz") into leg_ref_file
  script :
   headout="${params.output}_phas"
   """
   ${params.bin_bcftools} norm -m -any  $vcf -o tmp.vcf
    vcf2impute_legend_haps.pl -vcf tmp.vcf -leghap $headout -chr ${params.chr} #-start ${params.from_bp} -end ${params.to_bp}
   """
 }
}
process impute2{
  input :
     set file(haps_haps), file(haps_sampe) from file_haps_prephase
     file(maps) from genetic_map_ch_2
     file(haps_ref) from haps_ref_file
     file(leg_ref) from leg_ref_file
  publishDir "${params.output_dir}/impute2/intial", overwrite:true, mode:'copy'
  output :
    file("$headout*")
    set file("${headout}"), file(haps_sampe) into file_impute2
  script :
    size=params.to_bp!=""? " -int ${params.from_bp} ${params.to_bp} ": ""
    headout="${params.output}_impute2.gen"
    """
    ${params.bin_impute2}  ${size} \
      -known_haps_g ${haps_haps} \
      -h $haps_ref \
      -l $leg_ref \
      -m $maps \
      -Ne ${params.effective_size} \
      -buffer ${params.buffer_kb} \
      -o $headout
    """
}

process format_impute2{
  input :
     set file(genfile), file(sample) from file_impute2
  publishDir "${params.output_dir}/impute2/", overwrite:true, mode:'copy'
  output :
      file("$headout*")
      file("${headout}.vcf.gz") into (vcf_impute2gz, vcf_impute2gz_2)
  script :
    headout="${params.output}_impute2"
    chro=params.chro
    """
    awk -v chro=$chro '{\$2=chro":"\$3"_"\$4"_"\$5; print}' $genfile > $headout"_2.impute2"
    ${params.bin_bcftools} convert --gensample2vcf $headout"_2.impute2",$sample -o $headout".vcf"
    ${params.bin_bcftools} sort $headout".vcf"  -Oz -o  $headout".vcf.gz"
    """
}

process reag_imp2{
 input:
   file(impvcfgz) from vcf_impute2gz
 publishDir "${params.output_dir}/impute2/", overwrite:true, mode:'copy'
 output :
    file(headout)
 script :
   headout="${params.output}_impute2_reag.vcf"
   """
   ${params.bin_bcftools} norm -m +any  $impvcfgz -Oz -o $headout
   """
}




 genetic_map_ch_3=Channel.fromPath(params.genetic_map, checkIfExists:true)
 process postphase_shapeit{
  cpus params.cpus_other
  input :
      file(file_vcf)  from vcf_impute2gz_2
      file(genetic_map) from genetic_map_ch_3
  publishDir "${params.output_dir}/impute2/postphase", overwrite:true, mode:'copy'
  output :
      file("${fileout}*")
  script :
      chro=params.chr!=""? "--chrom=${params.chr}" : ""
      begin=params.to_bp!=""? "--bpStart=${params.from_bp}" : ""
      end=params.from_bp!=""? "--bpEnd=${params.to_bp}" : ""
      fileout="${params.output}_impute2"
      """
      ${params.bin_shapeit} \
       -V $file_vcf\
       --input-map $genetic_map \
       --input-thr ${params.thr_shapeit} \
       --output-max $fileout".haps" $fileout".sample" \
       --thread ${params.cpus_other} \
       --effective-size ${params.effective_size} \
       --output-log ${fileout}".log" \
       --force
        ${params.bin_shapeit} -convert \
        --input-haps  $fileout".haps" $fileout".sample" \
        --output-vcf $fileout".phased.vcf"
        ${params.bin_bcftools} norm -m +any  $fileout".phased.vcf" -Oz -o $fileout"_reagr.phased.vcf.gz"
      """
 }

