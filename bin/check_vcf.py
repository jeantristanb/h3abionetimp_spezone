#!/usr/bin/python
# coding=utf-8
import argparse
import re
import os
import sys
import glob
import gzip
def GetHeaderVcf(File):
    if readgz :
      VcfRead=gzip.open(File)
    else :
      VcfRead=open(File)
    Head=[]
    for line in VcfRead:
       if line[0]=="#" :
         Head.append(line.lower().replace("\n",""))
       else :
         VcfRead.close()
         return Head
    return Head

def checkrefalt(splline) :
  ref=splline[3]
  listalt=splline[4].split(',')
  if ('*' in listalt) or (ref in listalt) :
    return False
  return True


def parseArguments():
    parser = argparse.ArgumentParser(description='fill in missing bim values')
    parser.add_argument('--out',type=str,required=True)
    parser.add_argument('--vcf',type=str,required=True)
    args = parser.parse_args()
    return args


args = parseArguments()
typehead=args.vcf.split('.')[-1]
readgz=False
if typehead=='gz' or typehead=='.gzip':
  readgz=True
headervcf=GetHeaderVcf(args.vcf)
headervcf=headervcf[-1].split()
NCol=len(headervcf)



if readgz :
  readvcf=gzip.open(args.vcf)
else :
  readvcf=open(args.vcf)
writevcf=open(args.out, 'w')
for line in readvcf :
  if line[0]=="#" :
     writevcf.write(line)
  else :
   spll=line.split()
   if len(spll)!=NCol :
      print("\t".join(spll[0:5])+"not good size")
   elif checkrefalt(spll)==False :
        print(spll[1]+" "+spll[2]+": in ref "+spll[3]+" and alt "+spll[4])
   else :
       writevcf.write(line) 
