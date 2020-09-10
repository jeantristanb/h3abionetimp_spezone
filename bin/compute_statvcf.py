#!/usr/bin/python
# coding=utf-8
import argparse
import gzip
import re
import os
import sys
import glob





def get_info_pos(lines, postosave=None):
    def getgeno(x,l):
      if x[0]=='.' :
        return [None, None,0,0, None]
      x=x.split(':')[0]#re.split([':|/'])
      if '/' in x :
         x=x.split('/')
         phase=0
      elif '|' in x:
         x=x.split('|')
         phase=1
      return [l[int(x[0])], l[int(x[1])], int(x[0])+1, int(x[0])+1,phase]
    TmpSplit=lines.split()
    Chro=TmpSplit[0]
    Pos=TmpSplit[1]
    NomPos=Chro+"_"+Pos
    Ref=TmpSplit[3]
    ListeAlt=TmpSplit[4].split(',')
    TmpPos=[Ref]
    TmpPos+=ListeAlt
    if postosave :
     return (Chro, Pos, Ref, ListeAlt,[getgeno(TmpSplit[x], TmpPos)  for x in postosave ])
    else :
      return (Chro, Pos, Ref, ListeAlt,[getgeno(x, TmpPos)  for x in TmpSplit[9:]])

def get_header_vcf(File):
    VcfRead=open_file(File)
    Entete=[]
    for line in VcfRead:
       if line[0]=="#" :
         Entete.append(line.lower().replace("\n",""))
       else :
        VcfRead.close()
        return Entete
    VcfRead.close()
    return Entete

def open_file(File, Type='r'):
    typehead=File.split('.')[-1]
    readgz=False
    if typehead=='gz' or typehead=='.gzip':
      readgz=True
    try :
       if readgz :
         LireFich=gzip.open(File, Type)
       else :
         LireFich=open(File, Type)
    except IOError:
       sys.exit("File "+ File + " open in mode " + Type + "can't open File")
    return LireFich

def parseArguments():
    parser = argparse.ArgumentParser(description='fill in missing bim values')
    parser.add_argument('--out',type=str,required=True)
    parser.add_argument('--vcf',type=str,required=True)
    args = parser.parse_args()
    return args

def get_stat(loci, stat_ind, gethead=False):
 # (Chro, Pos, Ref, ListeAlt,[getgeno(x, TmpPos)  for x in TmpSplit[9:]])
 Head="chr\tpos\tref\talt\tnbalt\tdeletion\tnbinsertion\tnbrefhom\tnbalthom\tnbhet\tnbmissing"
 if gethead :
    return Head 
 ref=loci[2]
 nbalt=len(loci[3])
 delet=len(loci[2])
 insert=len([x for x in  loci[3] if len(x)>1])
 cmt_ind=0
 statsnp=[0,0,0,0] 
 for x in loci[4] :
  if x[0] :
    if x[0]!=x[1] :
       stat_ind[cmt_ind][2]+=1   
       statsnp[2]+=1
    elif x[0]==ref:
       stat_ind[cmt_ind][0]+=1   
       statsnp[0]+=1
    else :
       stat_ind[cmt_ind][1]+=1   
       statsnp[1]+=1
  else :
       stat_ind[cmt_ind][3]+=1   
       statsnp[3]+=1
  cmt_ind+=1 
 return loci[0]+"\t"+loci[1]+"\t"+ref+"\t"+",".join(loci[3])+"\t"+str(nbalt)+"\t"+str(delet)+"\t"+str(insert)+"\t"+"\t".join([str(x) for x in statsnp])



args = parseArguments()

headervcf=get_header_vcf(args.vcf)
NCol=len(headervcf[-1])
NameInd=headervcf[-1].split("\t")[9::]

readvcf=open_file(args.vcf)
statvcfind=[[0,0,0,0] for x in range(len(NameInd))]
writesnp=open(args.out+'_snp.txt', 'w')
writesnp.write(get_stat(None, None, True)+'\n')
for line in readvcf :
    if line[0]=="#" :
      continue 
    posinfo=get_info_pos(line) 
    writesnp.write(get_stat(posinfo, statvcfind)+'\n')
writesnp.close()
writeind=open(args.out+'_ind.txt', 'w')
Head="ID\tnbrefhom\tnbalthom\tnbhet\tnbmissing"
writeind.write(Head+'\n')
cmtind=0
for ind in statvcfind :
   writeind.write(NameInd[cmtind]+"\t"+"\t".join([str(x) for x in ind])+"\n")
writeind.close()
