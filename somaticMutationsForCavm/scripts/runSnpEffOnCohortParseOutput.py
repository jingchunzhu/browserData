#!/usr/bin/env python

import argparse
import os,sys
import re
import subprocess
import string
import gzip

os.sys.path.insert(0, os.path.dirname(__file__) )

def findSampleID (vcfFile): # this isRADIA specific
    fin =open(vcfFile,'r')
    while 1:
        line = fin.readline()
        if line[0]!="#":
            return ""
        if string.find(line,"##SAMPLE=<")!=-1 and ( string.find (line,"TUMOR")!=-1 or string.find(line,"Tumor") !=-1 or string.find(line,"tumor")!=-1 ) :# this isRADIA specific
            if string.find(line,"SampleTCGABarcode=")!=-1:
                id =string.split(string.split (line,"SampleTCGABarcode=")[1],",")[0]
                return id
        else:
            continue

def passingSomatic (vcfFile, directory, cavmId):
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
        elif string.find(line,"PASS")!=-1 and string.find (line,"SOM")!=-1  :# this isRADIA specific
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
    parser.add_argument("-passingSomatic", type=str, default="0",
                        help="1 for selecting passing somatics")
    args = parser.parse_args()

    fout = open(args.output,'w')
    fout.write("#"+string.join(["sample","chr","start","end","gene","reference","alt","effect","DNA_VAF","RNA_VAF","Amino_Acid_Change"],"\t")+"\n")
    fout.close()

    if args.id !="":
        cavmId = args.id
        vcfPathname = args.inputVcfDir
        if args.passingSomatic =="1":
            vcfPathname = passingSomatic(vcfPathname, "./", cavmId)

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
                if args.passingSomatic =="1":
                    vcfPathname = passingSomatic(vcfPathname, tmpDir, cavmId)
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
            cavmId= findSampleID (tmpDir+vcfFile)
            cmd= "cat "+ tmpDir+vcfFile +" | " +os.path.dirname(__file__)+ "/parseSnpEffVcf.py "+cavmId + " " + args.output
            subprocess.call(cmd, shell=True)

        os.system("rm -rf "+tmpDir)

if __name__ == '__main__':
        main()
    

