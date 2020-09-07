vcffile_ref=ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/working/20190425_NYGC_GATK/CCDG_13607_B01_GRM_WGS_2019-02-19_chr19.recalibrated_variants.vcf.gz
CrossMap=/home/jeantristan/Data/CrossMap/hg38ToHg19.over.chain.gz
chr=chr19
begin38=44899792
end38=44899826
around=100000
begin38r=`expr $begin38 \- $around`
end38r=`expr $end38 \+ $around`
fasta=~/Data/CrossMap/hg19.fa
input_col_ref=/spaces/jeantristan/Cassandra/1000GTest/Test_1000G_Cluster.xlsx
FileInd=/spaces/jeantristan/Cassandra/1000GTest/IndList_1000G.ind
genetic_map=/home/jeantristan/Data/Imputed2/1000GP_Phase3/genetic_map_chr19_combined_b37.txt
nextflow ../build_ref.nf --input_vcf_ref $vcffile_ref --maf 0.01 --chr $chr --from_bp $begin38r --to_bp $end38r -profile slurm --convert_file $CrossMap --fasta_file $fasta -resume --input_col_ref  $input_col_ref --keep $FileInd --genetic_map $genetic_map
#exit

## if file in bed 
GZF=All_20180423.vcf.gz ##
#nextflow run h3abionet/h3agwas/formatdata/plk_in_vcf_imp.nf --input_pat $FileBed2 --output_dir $chro"_"$begin38r"_"$end38r --file_ref_gzip $file_ref -resume --reffasta $fasta -profile slurm 
file_vcfi=reference_panel/infopos_add/refpanel_addpos.vcf
## for test we used subsample of snps using h3agwas position around 

begin37=45403041
end37=45403090
begin37r=`expr $begin37 \- $around`
end37r=`expr $end37 \+ $around`
awk -v begin=$begin37r -v end=$end37r '{if($1=="19" && $4>=(begin-100000) && $4<=(end+100000)){print "chr"$1"\t"$4"\t"$4}}' /dataE/AWIGenGWAS/plink/swt/input/swt.bim > pos_cheap

#nextflow  run ../imputed_data.nf --vcf $file_vcfi --maf 0.01 --chr $chr --from_bp $begin37r --to_bp $end37r  -profile slurm --bed pos_cheap --shapeit_ref reference_panel/shapeit/format/refpanel_phaseshapeit --refpanel_bref3 reference_panel/beagle/refpanel_phasebeagle.bref3
