#!/usr/bin/python
# coding=utf-8
import argparse
import re
import os
import sys
import glob


def GetHeaderVcf(File):
    VcfRead=open(File)
    Head=[]
    for line in VcfRead:
       if line[0]=="#" :
         Head.append(line.lower().replace("\n",""))
       else :
         VcfRead.close()
         return Head
    return Head

def parseArguments():
    parser = argparse.ArgumentParser(description='fill in missing bim values')
    parser.add_argument('--vcf',type=str,required=True)
    parser.add_argument('--out',type=str,required=True)
    args = parser.parse_args()
    return args

args = parseArguments()

headervcf=GetHeaderVcf(args.vcf)
headervcf=headervcf[-1].split()
Nind=len(headervcf)-9
readvcf=open(args.vcf)
list_pos=[]
list_typepos=[]
listinfoind=[]
for cmtind in range(Nind):
 listinfoind.append([[],[]]) 
rind=range(Nind)
nbpos=0
for line in readvcf :
   if line[0]!="#":
      splitl=line.split() 
      chro=splitl[0]
      list_pos.append(splitl[1])
      if splitl[4]==".":
        continue
      alt=splitl[4].split(",")
      if len(alt)==1 :
       list_typepos.append("S")
       error="?"
      else : 
       list_typepos.append("M")
       error="-1"
      for cmtind in rind:
        infoind=splitl[cmtind+9]
        if infoind[0]=='.' :
           geno=[error,error]
        else :
          geno=infoind.split(':')[0].split('/')    
        listinfoind[cmtind][0].append(geno[0]) 
        listinfoind[cmtind][1].append(geno[1]) 
      nbpos+=1
finalphase=str(Nind)+"\n"+str(nbpos)+"\nP "+" ".join(list_pos)+"\n"+" ".join(list_typepos)+"\n"
for cmtind in range(Nind) :
  finalphase+="#"+headervcf[cmtind+9]+"\n"
  finalphase+=" ".join(listinfoind[cmtind][0])+"\n"
  finalphase+=" ".join(listinfoind[cmtind][1])+"\n"

##
writePhase=open(args.out+'_PHASE.inp', 'w')
writePhase.write(finalphase)
writePhase.close()
