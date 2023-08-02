
/*
 * Authors       :
 *
 *
 *      jean-tristan Brandenburg
 *
 *  On behalf of the H3ABionet Consortium
 *  2015-2019
 *
 *
 * Description  : nextflow pipeline to phase data
 *
 *(C) University of the Witwatersrand, Johannesburg, 2020 
 *This is licensed under the MIT Licence. See the "LICENSE" file for details
 */

//---- General definitions --------------------------------------------------//

import java.nio.file.Paths;
import sun.nio.fs.UnixPath;
import java.security.MessageDigest;

allowed_params = ['file_vcf', 'ouptut_dir', 'output', "bin_PHASE"]

params.output_dir="out"
params.output="output"
params.file_vcf=""
params.bin_PHASE="PHASE"
params.bin_plink="plink"

if(params.file_vcf==""){
 error("\n\n------\nError in your config\nfile_vcf not defined\n\n---\n")
}

ch_filevcf=Channel.fromPath(params.file_vcf)

process formatdata{
  input:
   file(vcf) from ch_filevcf
  script :
    """ 
    plink --recode-fastphase

}
