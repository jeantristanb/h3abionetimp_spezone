## vcf initial
vcffile_ref=ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/working/20190425_NYGC_GATK/CCDG_13607_B01_GRM_WGS_2019-02-19_chr19.recalibrated_variants.vcf.gz
## to transform in hg19
CrossMap=/home/jeantristan/Data/CrossMap/hg38ToHg19.over.chain.gz
#chr and position to analyse
chr=chr19
begin38=44899792
end38=44899826
around=500000
begin38r=`expr $begin38 \- $around`
end38r=`expr $end38 \+ $around`

## fasta file correspondigng to transformation
fastahg19=~/Data/CrossMap/hg19.fa
input_col_ref=/spaces/jeantristan/Cassandra/1000GTest/Test_1000G_Cluster.xlsx
FileInd=/spaces/jeantristan/Cassandra/1000GTest/IndList_1000G.ind
genetic_map=/home/jeantristan/Data/Imputed2/1000GP_Phase3/genetic_map_chr19_combined_b37.txt

## build reference  for imputation 
newgenmap=chr19_genmap.genmap
if [ ! -f  $newgenmap ]
then
echo "chr position COMBINED_rate(cM/Mb) Genetic_Map(cM)" > $newgenmap
sed '1d' $genetic_map |awk -v chr=19 '{print chr" "$1" "$2" "$3" "$4}' >>  $newgenmap
fi

echo "nextflow ../build_ref.nf --input_vcf_ref $vcffile_ref --maf 0.01 --chr $chr --from_bp $begin38r --to_bp $end38r -profile slurm --convert_file $CrossMap --fasta_file $fastahg19 -resume --input_col_ref  $input_col_ref --keep $FileInd --genetic_map $newgenmap"
#exit

## if file in bed : transformation in vcf

## 
file_ref_gzip=All_20180423.vcf.gz ##
FileBedToImputed=""
#nextflow run h3abionet/h3agwas/formatdata/plk_in_vcf_imp.nf --input_pat $FileBedToImputed --output_dir "bed_in_vcf/" --file_ref_gzip $file_ref_gzip -resume --reffasta $fasta -profile slurm 
file_vcfi=reference_panel/vcfnorm/refpanel_sort_norm.vcf.gz
## for test we used subsample of snps using h3agwas position around 

##same position previous in hg19/37
around_imp=250000
begin37=45403041
end37=45403090
begin37r=`expr $begin37 \- $around_imp`
#end37r=`expr $end37 \+ $around_imp`
#begin37r=$begin37
end37r=`expr $end37 \+ $around_imp`
#end37r=$end37r
## to test we extracted some position of Array

#awk -v begin=$begin37r -v end=$end37r '{if($1=="19" && $4>=(begin-100000) && $4<=(end+100000)){print "chr"$1"\t"$4"\t"$4}}' /dataE/AWIGenGWAS/plink/swt/input/swt.bim > pos_cheap

## imputation with beagle
#nextflow  run ../imputed_data_beagle.nf --vcf $file_vcfi --maf 0.01 --chr $chr --from_bp $begin37r --to_bp $end37r  -profile slurm --bed pos_cheap --vcf_ref reference_panel/beagle/refpanel_phasebeagle.vcf.gz  
nextflow  run ../imputed_data_impute2.nf --vcf $file_vcfi --maf 0.01 --chr $chr --from_bp $begin37r --to_bp $end37r  -profile slurm --bed pos_cheap --vcf_ref reference_panel/shapeit/vcf/refpanel_phaseshapeit.phased.vcf --genetic_map $genetic_map -resume

