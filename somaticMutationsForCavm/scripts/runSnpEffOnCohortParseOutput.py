#!/usr/bin/env python

import argparse
import os,sys
import re
import subprocess
import string

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

def passingSomatic (vcfFile, cavmId):
    ## really stupid thing to get rid of vcf headers
    output = ".tmp_"+cavmId
    fin =open(vcfFile,'r')
    fout= open(output,'w')
    while 1:
        line = fin.readline()
        if line=="":
            break
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
            vcfPathname = passingSomatic(vcfPathname, cavmId)

        cmd = "export PATH="+ os.path.dirname(__file__)+"/:$PATH; runSnpEffAgainstRefSeq.bash "+vcfPathname +" | "+ os.path.dirname(__file__)+ "/parseSnpEffVcf.py "+cavmId + " " + args.output
        subprocess.call(cmd, shell=True)

        if vcfPathname != args.inputVcfDir:
            os.system("rm "+vcfPathname)

    else:
        for vcfFile in os.listdir(args.inputVcfDir):
            if re.search("\.vcf$", vcfFile):
                #file size >0
                if args.inputVcfDir[-1]=="/":
                    if os.stat(args.inputVcfDir+vcfFile).st_size ==0:
                        continue
                else:
                    if os.stat(args.inputVcfDir+"/"+vcfFile).st_size ==0:
                        continue
                if args.id != "":
                    cavmId = args.id
                else:
                    cavmId= findSampleID (args.inputVcfDir+"/"+vcfFile)
                
                if cavmId =="":
                    print "no id specified in vcf file hearder or through command line\n"
                    sys.exit()

                vcfPathname = args.inputVcfDir + "/" + vcfFile

                #passing somatic
                if args.passingSomatic =="1":
                    vcfPathname = passingSomatic(vcfPathname, cavmId)
                
                print vcfPathname
                cmd = "export PATH="+ os.path.dirname(__file__)+"/:$PATH; runSnpEffAgainstRefSeq.bash "+vcfPathname +" | "+ os.path.dirname(__file__)+ "/parseSnpEffVcf.py "+cavmId + " " + args.output
                subprocess.call(cmd, shell=True)
        
                if vcfPathname != args.inputVcfDir + "/" + vcfFile: 
                    os.system("rm "+vcfPathname)

if __name__ == '__main__':
        main()
    

