import string, os,sys
import parseSnpEffVcf

def getSeq_hg19(chr,start,end):
    hg19_fasta= "/inside/depot4/CCI/human_g1k_v37.fasta"  ######### hard coded for hg19
    os.system("samtools faidx " + hg19_fasta +" "+chr+":"+str(start)+"-"+str(end) +" | tail -n 1 > .seq")
    fin = open(".seq",'r')
    seq = string.strip(fin.read())
    fin.close()
    return seq


#parse SNPeff output and remake xena file
def parseSNPeffToXenaOutput(filein, fileout, xenaDic):
    fout = open(fileout,'w')
    fin = open(filein,'r')
    myVcf = parseSnpEffVcf.vcf(fin)
    fin.close()

    for row in myVcf.read():
        id = string.join([row.ID,row.chr[3:],str(row.start)],"_")
        sample, chr, start, end, oldgene, ref,alt, var_class, DNA_AF, RNA_AF, AA_change = xenaDic[id]
        if len(row.effectPerGene)!=0:
            for gene in row.effectPerGene.keys():
                AA_Change = row.effectPerGene[gene]["Amino_Acid_Change"]
                if string.find(row.effectPerGene[gene]["effect"],"frameshift")!=-1:
                    AA_Change =""
                fout.write(string.join([sample, row.chr, str(row.start),
                                        str(row.end), gene, ref, alt, 
                                        row.effectPerGene[gene]["effect"], DNA_AF, RNA_AF, AA_Change],"\t")+"\n")
        else:
            fout.write(string.join([sample, row.chr, str(row.start),
                                    str(row.end), "", ref, alt, 
                                    "", DNA_AF, RNA_AF, ""],"\t")+"\n")
    fout.close()

if len(sys.argv[:])!= 3:
    print "python xenaReAnnotateHg19.py xenaIn xenaNewOut\n"
    sys.exit()

fin = open(sys.argv[1],'r')
vcfFile = ".tmp"
snpeffFile = ".tmp.snpeff"
output =sys.argv[2]

#build vcf file
dic = {}
fvcf = open(vcfFile,'w')
for line in fin.readlines():
    if line[0]=="#":
        continue
    line =line[:-1]
    data =string.split(line,'\t')
    sample, chr, start, end, gene, ref,alt, var_class, DNA_AF, RNA_AF, AA_change =data
    id = string.join([sample,chr,start],"_")

    if id in dic:
        print "problem with ", id
        continue
    if ref in ["-","."]:
        ref = getSeq_hg19(chr, int(start), int(start))
        alt = alt+ref
        end =start
    elif alt in ["-","."]:
        nextBase = getSeq_hg19(chr, int(start)+1,int(start)+1) 
        alt=nextBase
        ref = ref+ nextBase
        end =str (int(start) + len(ref) -1)

    dic[id]=[sample, chr, start, end, gene, ref,alt, var_class, DNA_AF, RNA_AF, AA_change]
    fvcf.write(string.join([chr,start,sample, ref,alt],'\t')+'\n')
fvcf.close()

parseSNPeffToXenaOutput (snpeffFile, output, dic)
sys.exit()
#re-run SNPeff
dir ="/inside/home/jzhu/scripts/vcfXenaData/browserDataMelisssa/somaticMutationsForCavm/"
os.system("java -jar -Xmx2480m "+dir + "/snpEff/snpEff.jar hg19 -c "+ dir + "snpEff/snpEff.config -quiet -sequenceOntology -no-downstream -no-upstream -no-intergenic -no-intron -v -noStats .tmp > "+ snpeffFile )

