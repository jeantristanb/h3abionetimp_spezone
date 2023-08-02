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

allowed_params = ["input_plk_toimp", "input_vcf_toimp","inputs_ref"]

params.each { parm ->
  if (! allowed_params.contains(parm.key)) {
    println "\nUnknown parameter : Check parameter <$parm>\n";
  }
}


def params_help = new LinkedHashMap(helps)


params.queue      = 'batch'
params.work_dir   = "$HOME/h3agwas"
params.output_dir = "${params.work_dir}/output"

params.input_plk_toimp=""
params.input_vcf_toimp=""
params.inputs_ref=""

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


if(params.input_plk_toimp=="" && params.input_vcf_toimp==""){
  println "This workflow requires input_plk_toimp or input_vcf_toimp"
  exit 1
}
if(params.inputs_ref==""){
  println "This workflow requires inputs_ref"
  exit 1

}


if(params.input_plk_toimp!="" && params.input_vcf_toimp!=""){
  println "This workflow requires input_plk_toimp or input_vcf_toimp but not both : \n\tinput_plk_toimp : ${params.input_plk_toimp} \n\t input_vcf_toimp : ${params.input_vcf_toimp}"
  exit 1
}


/*Step : build with */
if(params.input_plk_toimp!=""){


}
