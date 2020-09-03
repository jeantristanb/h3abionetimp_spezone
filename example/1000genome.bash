vcffile=ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/working/20190425_NYGC_GATK/CCDG_13607_B01_GRM_WGS_2019-02-19_chr19.recalibrated_variants.vcf.gz
CrossMap=/home/jeantristan/Data/CrossMap/hg38ToHg19.over.chain.gz
chr=chr19
begin=44899792
end=44899826
around=100000
beginr=`expr $begin \- $around`
endr=`expr $end \+ $around`
fasta=~/Data/CrossMap/hg19.fa
nextflow /home/jeantristan/Travail/git/PhaseImp/h3abionetimp_spezone/build_ref.nf --input_vcf_ref $vcffile --maf 0.01 --chr $chr --from_bp $beginr --to_bp $endr -profile slurm --convert_file $CrossMap --fasta_file $fasta -resume
