#!/usr/bin/env Rscript
library("readxl")
getallele<-function(rs, chro,  position,Id,InfoVcf,a1,a2, filecorresout){
Id<-as.character(Id);a1<-as.character(a1);a2<-as.character(a2);rs=as.character(rs);chro=as.character(chro);position=as.character(position)
NewList<-c('A','T','G','C', 'AA','AT', 'AC', 'AG',  'TA','TT', 'TC', 'TG', 'CA','CT', 'CC', 'CG', 'GG','GT', 'GC', 'GG')
QUAL="."
Filter="PASS"
INFO="."
FORMAT="GT"
a1<-toupper(a1)
a2<-toupper(a2)
a1[is.na(a1)]<-'.'
a2[is.na(a2)]<-'.'
tba1a2<-unique(toupper(c(a1,a2)))
perc<-sapply(strsplit(tba1a2[tba1a2!='.'], split=''), function(x)length(which(x %in% c('A','T', 'C', 'G')))/length(x))
tba1a2nomiss<-tba1a2[tba1a2!='.']
nballele<-length(tba1a2nomiss)
if(all(perc==1)){
corr_allele<-data.frame(InitialAll=c('.',tba1a2nomiss), UseAll=c('.',tba1a2nomiss),NumAll=c('.',0:(nballele-1)), stringsAsFactors=F)
}else{
corr_allele<-data.frame(InitialAll=c('.',tba1a2nomiss), UseAll=c('.',NewList[1:nballele]),NumAll=c('.',0:(nballele-1)), stringsAsFactors=F)
}
write.csv(corr_allele, file=paste(filecorresout,'_',rs,'_corrallele.csv',sep=''),row.names=F)
Ref<-corr_allele[2,2]
Alt<-paste(corr_allele[3:nrow(corr_allele),2],collapse=',')
Allele<-merge(data.frame(Id=Id ,A1=a1, A2=a2, stringsAsFactors=F), InfoVcf, by=1,all.y=T)
Allele$A1[is.na(Allele$A1)]<-'.'
Allele$A2[is.na(Allele$A2)]<-'.'
tmpcorr<-corr_allele[,c('InitialAll', 'NumAll')]
names(tmpcorr)<-c('A1', 'NumA1')
Allele1<-merge(Allele, tmpcorr, by='A1',all.x=T)
names(tmpcorr)<-c('A2', 'NumA2')
Allele2<-merge(Allele1, tmpcorr, by='A2',all.x=T)
Allele2$GenoVcf<-paste(Allele2$NumA1,Allele2$NumA2 ,sep='/')
Allele2<-Allele2[order(Allele2$Num),]
PosVcf<-paste(c(chro, position, rs, Ref, Alt, QUAL, Filter, INFO,FORMAT, Allele2$GenoVcf), collapse="\t")
return(PosVcf)
}

args = commandArgs(trailingOnly=TRUE)
#FileTomm='vcf_wgs_APOE_20200601.xlsx'
if(length(args)==0){
FileTom="/home/jeantristan/TOMM40_Dist.xlsx"
FileVcf="/spaces/jeantristan/Cassandra/MergeWithToom40/2_MergeDataI/19_45353041_45453090.recode.vcf"
FileOut="test"
}else{
FileTom=args[1]
FileVcf=args[2]
FileOut=args[3]
}

if(length(grep('.xls$', FileTom))==1 | length(grep('.xlsx$', FileTom))==1){
DataI<-as.data.frame(read_excel(FileTom))
}else if (grep('.csv$', FileTom)==1){
DataI<-read.csv(FileTom)
}else{
DataI<-read.table(FileTom, header=T)
}

ListRS=names(DataI)[-1]
Tmp<-grep("^#",readLines(FileVcf, 5000),value=T)
HeadName<-Tmp[length(Tmp)]
ListIndTmp<-strsplit(HeadName,split='\t')[[1]]
ListIndVcf<-ListIndTmp[10:length(ListIndTmp)]
PosIndVcf<-data.frame(NameInd=ListIndVcf, Num=1:length(ListIndVcf))
if(length(which(DataI[,1] %in% ListIndVcf))==0){
cat('error\ntransform_file_in_vcf.r : no name in common with vcf\n')
cat('vcf list :',paste(ListIndVcf,collapse=','), '\n')
cat('DataI :',paste(as.character(DataI[,1]),collapse=','), '\n')
q(status = 2)
}

##
baliseind<-DataI[,1] %in% ListIndVcf
MatRs<-cbind(as.data.frame(matrix(unlist(strsplit(names(DataI)[-1],split='_')),ncol=4,byrow=T)), Head=names(DataI)[-1])
ListeAllFormVcf<-c(HeadName)
for(rs in unique(MatRs[,1])){
MatRsTmp<-MatRs[MatRs[,1]==rs,]
if(nrow(MatRsTmp)!=2){
print(MatRsTmp)
cat('error ', rs, ' more than one observation\n')
q(status=2)
}
HeadA1<-as.character(MatRsTmp[1,5])
HeadA2<-as.character(MatRsTmp[2,5])
#getallele<-function(rs, chro,  position,Id,InfoVcf,a1,a2, filecorresout){
chro<-unique(MatRsTmp[,2]);position<-unique(MatRsTmp[,3])
ListeAllFormVcf<-c( ListeAllFormVcf,getallele(rs,chro,  position,DataI[baliseind, 1],PosIndVcf,DataI[baliseind, HeadA1],DataI[baliseind, HeadA2],FileOut))
}
writeLines(ListeAllFormVcf, con=paste(FileOut,'.vcf',sep=''))
