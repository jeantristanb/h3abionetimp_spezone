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

// Checks if the file exists
checker = { fn ->
   if (fn.exists())
       return fn;
    else
       error("\n\n-----------------\nFile $fn does not exist\n\n---\n")
}

import java.nio.file.Paths



def helps = [ 'help' : 'help' ]

allowed_params = ["input_vcf_ref", "input_col_ref", "keep", "chr", "from_bp", "to_bp", "extract", "keep", "output_dir", "maf", "output", "big_time", "convert_file", "fasta_file"]
infosys=["scripts"]
infosoft=["bin_vcftools", "bin_vcftools", "big_time", "bin_tabix", "bin_crossmap", "bin_bcftools", "bin_beagle", "bin_shapeit", "bin_bref3"]
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
params.convert_file=""
/**/
//params.bp_flanking=1

/**/
params.thr_shapeit=0.95
/*software information */

params.bin_beagle="beagle"
params.bin_tabix="tabix"
params.bin_vcftools="vcftools"
params.bin_bcftools="bcftools"
params.bin_shapeit="shapeit"
params.bin_bref3="bref3"
params.bin_eagle="eagle"

params.memory_vcftools="10GB"
params.memory_tabix="10GB"
params.cpus_other=4
params.bin_crossmap="~/.local/bin/CrossMap_beta.py"
params.exclude_bed=""
params.bed=""
params.genetic_map=""
params.genetic_map_beagle=""
params.max_missing=0.99

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


vcfrefisgz=file_vcf_tofilter.getExtension()=='gz'
}else{
      file_vcf_tofilter=Channel.fromPath(params.input_vcf_ref)
      vcfrefisgz=file(params.input_vcf_ref).getExtension()=='gz'
}

//case where need to do some filtering or added information
if(params.input_col_ref!="" || params.keep!="" || params.chr!="" || (params.chr!="" && params.to_bp!="" && params.from_bp!="") || params.keep!="" || params.maf!="" || params.convert_file!="" || params.exclude_bed!=""|| params.max_missing<1){

  if(params.keep!="" || params.chr!="" || (params.chr!="" && params.to_bp!="" && params.from_bp!="") || params.exclude_bed!="" || params.bed!="" || params.max_missing<1){
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
      output :
         file("${out_file}.recode.vcf") into file_vcf_filter_1
      script :
         gvcf= file_vcf.getExtension()  ? "--gzvcf" : "--vcf" 
         chro=params.chr!=""? "--chr ${params.chr}" : ""
         end=params.to_bp!=""? "--to-bp  ${params.to_bp}" : "" 
         begin=params.from_bp!=""? "--from-bp  ${params.from_bp}" : "" 
         maf=params.maf!="" ? "--maf ${params.maf} " : ""
         keep=params.keep!="" ? " --keep ${file_ind} " : ""
         excl_bed=params.exclude_bed!="" ? " --exclude_bed $exclude_bed " : "" 
         bed=params.bed!="" ? " --bed $bedf " : ""
         out_file="${params.output}_filt1"
         maxmissing=params.max_missing<1 ? " --max-missing ${params.max_missing} " : ""
         """
         ${params.bin_vcftools} $gvcf $file_vcf $chro $end $begin $maf $keep --out ${out_file} --recode --recode-INFO-all $excl_bed $bed $maxmissing
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
    vcfrefisgz=0

  }else{
    file_vcf_filter_2 = file_vcf_filter_1
  }
  if(params.input_col_ref!=""){
     inputcolref=Channel.fromPath(params.input_col_ref)
     process AddPosition {
        input:
          file(filevcfi) from file_vcf_filter_2
          file(colreffile) from inputcolref  
         output :
           file("${fileposvcf}*")
           file("$filevcffinal") into file_vcf_filter_final 
        script :
          fileposvcf="${params.output}_pos"
          filevcffinal="${params.output}_addpos.vcf"
          readvcf=(vcfrefisgz)? "zcat " : "cat"
          """ 
          $readvcf $filevcfi > $filevcffinal
          transform_file_in_vcf.r $colreffile $filevcfi $fileposvcf
          grep -v '#' $fileposvcf".vcf" >> $filevcffinal 
          """
     }   
     vcfrefisgz=0
  }else{
    file_vcf_filter_final=file_vcf_filter_2
  }

}else{
  file_vcf_filter_final=file_vcf_tofilter
}

process IndexFile{
  input :
     file(filevcf) from file_vcf_filter_final
  publishDir "${params.output_dir}/vcffilt/", overwrite:true, mode:'copy'
  output :
     set file("${fileout}.gz"), file("${fileout}.gz.csi") into (file_vcf_filter_index,file_vcf_filter_index2)
  script :
      fileout="${params.output}_sort.vcf"
      logfilt="${params.output}_filter.log"
      readvcf=(vcfrefisgz)? "zcat " : "cat"
      """
      check_vcf.py --vcf $filevcf --out vcftmp.vcf > $logfilt
      ${params.bin_bcftools} sort vcftmp.vcf| bgzip -c > $fileout".gz"
      ${params.bin_bcftools} index $fileout".gz"
      """
}


process  normvcf{
  input :
     set file(filevcf), file(filevcfidx) from file_vcf_filter_index2
  publishDir "${params.output_dir}/vcfnorm/", overwrite:true, mode:'copy'
  output :
     set file("${fileout}"), file("${fileout}.csi") into (file_vcf_filter_index_norm , file_vcf_filter_index_norm_2, file_vcf_filter_index_norm_3)
  script :
      fileout="${params.output}_sort_norm.vcf.gz" 
      """
      ${params.bin_bcftools} norm -m -any $filevcf -Oz -o $fileout
      ${params.bin_bcftools} index $fileout
      """
}

if(params.genetic_map_beagle!=""){
      gm_ch=Channel.fromPath(params.genetic_map)
     }else{
      gm_ch=file('nogenetmap')
}

process phasebeagle{
 cpus params.cpus_other
 input :
   file(genet_map) from gm_ch
   set file(filevcf), file(filevcfidx) from file_vcf_filter_index
 publishDir "${params.output_dir}/beagle/", overwrite:true, mode:'copy'
 output :
   file("${headout}*")
   set file("${headout}.vcf.gz"), file("${headout}.vcf.gz.csi") into file_vcf_beaglephase
 script :
   map=(params.genetic_map_beagle!="")?"map=$genet_map" : ""
   headout="${params.output}_phasebeagle"
   """
   ${params.bin_beagle} gt=$filevcf nthreads=${params.cpus_other} out=$headout  $map ${params.otheropt_beagle}
   ${params.bin_bcftools} index $headout".vcf.gz"
   """
}

process convertbeagle{
 input :
   set file(filevcf), file(filevcfidx) from file_vcf_beaglephase
 publishDir "${params.output_dir}/beagle/", overwrite:true, mode:'copy'
 output :
   file("${headout}")
 script :
   headout="${params.output}_phasebeagle.bref3"
   """
   ${params.bin_bref3}  $filevcf > $headout
   """
}




gm_ch2=Channel.fromPath(params.genetic_map)
process shapeit_phase{
 cpus params.cpus_other
 input :
   file(genet_map) from gm_ch2
   set file(filevcf), file(filevcfidx) from file_vcf_filter_index_norm
 publishDir "${params.output_dir}/shapeit/init/", overwrite:true, mode:'copy'
 output :
   file("${headout}*") 
   set file("${headout}_i.haps"), file("${headout}_i.sample") into (shapeit_init_1, shapeit_init_2)
 script :
   headout="${params.output}_phaseshapeit"
   """
    awk '{print \$2" "\$3" "\$4}' $genet_map > tmpmap.txt
   ${params.bin_shapeit} -V $filevcf --input-map tmpmap.txt --input-thr ${params.thr_shapeit} --output-max $headout"_i.haps" $headout"_i.sample" --thread ${params.cpus_other}  --force

   """
}
process shapeit_convert_haps{
 input :
  set file(haps), file(sample) from shapeit_init_1
 publishDir "${params.output_dir}/shapeit/format/", overwrite:true, mode:'copy'
 output :
  file("$headout*")
 script : 
   headout="${params.output}_phaseshapeit"
   """
   ${params.bin_shapeit} -convert --input-haps  $haps $sample --output-ref  $headout".haps" $headout".leg" $headout".sample"
   """
}


process shapeit_convert_vcf{
 input :
  set file(haps), file(sample) from shapeit_init_2
 publishDir "${params.output_dir}/shapeit/vcf/", overwrite:true, mode:'copy'
 output :
  file("$headout*")
 script :
   headout="${params.output}_phaseshapeit"
   """
   ${params.bin_shapeit} -convert \
        --input-haps  $haps $sample \
        --output-vcf $headout".phased.vcf"
   """
}


gm_ch3=Channel.fromPath(params.genetic_map)
process eagle_phase{
   input :
     set file(filevcf), file(filevcfidx) from file_vcf_filter_index_norm_3
     file(map) from gm_ch3
    publishDir "${params.output_dir}/eagle/vcf/", overwrite:true, mode:'copy'
    output :
      file("${file_out}*")
    script :
      chro=params.chr!=""? "--chrom=${params.chr}" : ""
      begin=params.to_bp!=""? "--bpStart=${params.to_bp}" : ""
      end=params.from_bp!=""? "--bpEnd=${params.from_bp}" : ""
      file_out="${params.output}_prephase_eagle"
      """
        ${params.bin_eagle}\
                --vcf=${filevcf} \
                --geneticMapFile=${map} \
                --vcfOutFormat=z \
                $chro \
                $begin \
                $end \
                --maxMissingPerSnp=${params.max_missing} \
                --maxMissingPerIndiv=${params.max_missing} \
                --outPrefix=${file_out} 2>&1 | tee ${file_out}.log
      """
}

