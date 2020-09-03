wget http://stephenslab.uchicago.edu/assets/software/phase/phasecode/phase.2.1.1.linux.tar.gz
tar -xzf phase.2.1.1.linux.tar.gz
mv phase.2.1.1.linux/PHASE .
chmod +x PHASE
rm -rf phase.2.1.1.linux phase.2.1.1.linux.tar.gz
#git clone https://github.com/lstevison/vcf-conversion-tools
#cp vcf-conversion-tools/vcf2fastPHASE.pl .
#rm -rf vcf-conversion-tools
#chmdo +x vcf2fastPHASE.pl
wget https://mathgen.stats.ox.ac.uk/genetics_software/shapeit/shapeit.v2.r904.glibcv2.12.linux.tar.gz 
tar -xzf shapeit.v2.r904.glibcv2.12.linux.tar.gz
cp shapeit.v2.904.2.6.32-696.18.7.el6.x86_64/bin/shapeit .
rm -rf shapeit.v2.904.2.6.32-696.18.7.el6.x86_64/ shapeit.v2.r904.glibcv2.12.linux.tar.gz
