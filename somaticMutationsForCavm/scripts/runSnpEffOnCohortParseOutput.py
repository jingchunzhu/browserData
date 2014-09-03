#!/usr/bin/env python

import argparse
import os,sys
import re
import subprocess
import string
import gzip

os.sys.path.insert(0, os.path.dirname(__file__) )

def findSampleID (vcfFile): # this isRADIA specific
    if re.search("\.vcf.gz$", vcfFile):
        fin = gzip.open(vcfFile, 'rb')
    elif re.search("\.vcf$", vcfFile):
        fin =open(vcfFile,'r')
    while 1:
        line = fin.readline()
        if line[0]!="#":
            return ""
        #SAMPLE=<ID=DNA_NORMAL,...>
        if string.find(line,"##SAMPLE=<")!=-1 :
            found =0
            pairs = string.split(string.split(line,"##SAMPLE=<")[1][:-1],",")
            for pair in pairs:
                key,value=string.split(pair,"=")[0:2]
                #RADIA
                if key=="ID" and value=="DNA_TUMOR":
                    found =1
                    break
            if not found:
                continue

            for pair in pairs:
                key,value=string.split(pair,"=")[0:2]
                #TCGA
                if key=="SampleTCGABarcode":
                    return value
        else:
            continue

def passingSomatic (vcfFile, directory, cavmId, code):
    output = directory+".tmp_"+cavmId
    if re.search("\.vcf.gz$", vcfFile):
        fin = gzip.open(vcfFile, 'rb')
    elif re.search("\.vcf$", vcfFile):
        fin =open(vcfFile,'r')

    fout= open(output,'w')
    while 1:
        line = fin.readline()
        if line=="":
            break
        if line[0]=="#":
            fout.write(line)
        elif code==1 and string.find(line,"PASS")!=-1 and string.find (line,"SOM") !=-1 and string.find(line,"SS=2")!=-1  :# passing somatic
            fout.write(line)
        elif code==2 and string.find(line,"PASS")!=-1 and string.find (line,"SOM")!=-1 and  string.find(line,"SS=2")  and string.find (line,"VT=SNP")!=-1 :# passing somatic SNP
            fout.write(line)
        else:
            continue
    fout.close()
    return output

    
def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("inputVcfDir", type=str,
                        help="Directory of input VCFs, with wildcards | a single VCF file")
    parser.add_argument("output", type=str,
                        help="Xena Output file")
    parser.add_argument("-id", type=str, default="",
                        help="CAVM ID")
    parser.add_argument("-passingSomatic", type=int, default=0,
                        help="1 for selecting passing somatics, 2 for selecting passing somatic SNPs")
    args = parser.parse_args()

    fout = open(args.output,'w')
    fout.write("#"+string.join(["sample","chr","start","end","reference","alt","gene","effect","DNA_VAF","RNA_VAF","Amino_Acid_Change"],"\t")+"\n")
    fout.close()

    if args.id !="":
        cavmId = args.id
        vcfPathname = args.inputVcfDir
        if args.passingSomatic >0 :
            vcfPathname = passingSomatic(vcfPathname, "./", cavmId, args.passingSomatic)

        cmd = "export PATH="+ os.path.dirname(__file__)+"/:$PATH; runSnpEffAgainstRefSeq.bash "+vcfPathname +" | "+ os.path.dirname(__file__)+ "/parseSnpEffVcf.py "+cavmId + " " + args.output
        subprocess.call(cmd, shell=True)

        if vcfPathname != args.inputVcfDir:
            os.system("rm "+vcfPathname)

    else:
        fList=open(".fileList",'w')
        tmpDir="new/"
        if os.path.exists(tmpDir):
            os.system("rm -rf "+ tmpDir)
        os.system("mkdir "+tmpDir)


        for vcfFile in os.listdir(args.inputVcfDir):
            if re.search("\.vcf$", vcfFile) or re.search("\.vcf.gz$", vcfFile):
                #file size >0
                if args.inputVcfDir[-1]=="/":
                    if os.stat(args.inputVcfDir+vcfFile).st_size ==0:
                        continue
                else:
                    if os.stat(args.inputVcfDir+"/"+vcfFile).st_size ==0:
                        continue
 
                vcfPathname = args.inputVcfDir + "/" + vcfFile
                cavmId= findSampleID (args.inputVcfDir + "/" + vcfFile)

                #passing somatic
                if args.passingSomatic >0 :
                    vcfPathname = passingSomatic(vcfPathname, tmpDir, cavmId, args.passingSomatic)
                    fList.write(vcfPathname+"\n")
                elif re.search("\.vcf$", vcfFile):
                    os.system("cp " +  args.inputVcfDir + "/" + vcfFile + " "+ tmpDir)
                    fList.write(tmpDir+vcfFile+"\n")
                elif re.search("\.vcf.gz$", vcfFile):
                    os.system("cp " +  args.inputVcfDir + "/" + vcfFile + " "+ tmpDir + "| gunzip "+ tmpDir+vcfFile )
                    fList.write(tmpDir+vcfFile[:-3]+"\n")
            
        fList.close()
        cmd = "export PATH="+ os.path.dirname(__file__)+"/:$PATH; runSnpEffAgainstRefSeqFileList.bash .fileList"
        subprocess.call(cmd, shell=True)

        for vcfFile in os.listdir(tmpDir):
            if re.search("\.eff.vcf$", vcfFile) :
                cavmId= findSampleID (tmpDir+vcfFile)
                cmd= "cat "+ tmpDir+vcfFile +" | " +os.path.dirname(__file__)+ "/parseSnpEffVcf.py "+cavmId + " " + args.output
                subprocess.call(cmd, shell=True)

        os.system("rm -rf "+tmpDir)

if __name__ == '__main__':
        main()
    

