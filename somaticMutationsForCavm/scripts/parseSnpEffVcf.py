#!/usr/bin/env python

import argparse
import re
import sys, os
import string
import math

impact = {"MODIFIER":0, "LOW":1, "MODERATE":2, "HIGH":3}

class vcfRow(object):
    """This object contains an individual row of VCF output"""
    def __init__(self, row, columnLabels, EFFECT=1):
        tokens = row[:-1].split("\t")
        #chromosome notation
        chrom = tokens[0]
        if string.upper(chrom[0:3])== "CHR" and chrom[0:3]!="chr":
            chrom="chr"+chrom[3:]
        elif string.upper(chrom[0:2])== "CH" and string.upper(chrom[2])!="R":
            chrom= "chr"+chrom[2:]
        elif chrom in ["23","X","x"]:
            chrom="chrX"
        elif chrom in ["24","Y","y"]:
            chrom="chrY"
        elif chrom in ["25","M","m"]:
            chrom="chrM"
        else:
            chrom="chr"+chrom

        if chrom == "chr23":
            chrom="chrX"
        if chrom == "chr24":
            chrom="chrY"

        self.chr = chrom
        self.start = int(tokens[1])
        self.ID = tokens[2]
        self.reference = tokens[3]
        self.alt = tokens[4]
        self.end = self.start + len(self.reference) - 1
        self.DNA_AF=""
        self.RNA_AF=""
        self.NORMAL_AF=""

        effectsString = ""
        for thisSubToken in tokens[7].split(";"):
            if re.search("^EFF=", thisSubToken):
                effectsString = thisSubToken

        GT_code=None

        if len(tokens)<=8:
            pass
        else:
            # this could be all specific to RADIA output
            format =tokens[8]
            DNA_NORMAL = tokens[9]
            DNA_TUMOR =tokens[10]
            if len(tokens)>11:
                RNA_TUMOR = tokens[11]
            else:
                RNA_TUMOR=''

            #get the ALT base code (in VCF: GT=0,1,2?, 1 and 2 are acceptable)
            GT_code = self._findGTCode(self.chr, format, DNA_TUMOR, RNA_TUMOR, self.start)

            if GT_code !=None:
                self.alt = string.split(tokens[4],",")[GT_code-1]
                # AF allel frequency # this is all specific to RADIA output
                ID="AD"
                val =self._parse_TUMOR_ALT_ID (ID,format,DNA_TUMOR, RNA_TUMOR, GT_code)
                if val != None:
                    self.DNA_AD, self.RNA_AD= val
                    ID="DP"
                    val = self._parse_TUMOR_SINGLE_ID(ID,format,DNA_TUMOR, RNA_TUMOR)
                    if val != None:
                        self.DNA_DP, self.RNA_DP= val
                        try:
                            self.DNA_AF = str(float(self.DNA_AD) /float(self.DNA_DP))
                        except:
                            self.DNA_AF = "NA"
                        try:
                            self.RNA_AF = str(float(self.RNA_AD) /float(self.RNA_DP))
                        except:
                            self.RNA_AF = "NA"
                    else:
                        self.DNA_AF, self.RNA_AF= "NA","NA"
                else: # returned None
                    #print "AF error", row
                    self.DNA_AF, self.RNA_AF= "NA","NA"

                #get info on normal sample
                ID = "AD"
                val = self._parse_NORMAL_ALT_ID (ID,format,DNA_NORMAL, GT_code)
                if val != None:
                    if val =="NA":
                        val =0
                    self.NORMAL_AD = val
                    ID ="DP"
                    val = self._parse_NORMAL_SINGLE_ID(ID,format,DNA_NORMAL)
                    if val != None:
                        self.NORMAL_DP = val
                        try:
                            self.NORMAL_AF = str(float(self.NORMAL_AD) /float(self.NORMAL_DP))
                        except:
                            self.NORMAL_AF = "NA"
                    else:
                        self.NORMAL_AF = "NA"
                else:
                    self.NORMAL_AF = "NA"
            else:
                #print "GT error", row
                self.alt = "NA"
                self.DNA_AF, self.RNA_AF= "NA","NA"


        if EFFECT:
            self.effectPerGene = self._parseEffectsPerGene(effectsString, columnLabels, GT_code)
        return

    def get_DNAVAF(self):
        return self.DNA_AF

    def get_RNAVAF(self):
        return self.RNA_AF

    def get_NORMALVAF(self):
        return self.NORMAL_AF

    def _findGTCode(self, chrom, format, DNA_TUMOR, RNA_TUMOR, start):
        pos=-1
        data= string.split(format,":")
        for i in range(0,len(data)):
            if data[i]== "GT":
                pos=i
                break
        if pos ==-1:
            return None

        DNA_AD=-1
        RNA_AD=-1
        DNA_GT_code =-1
        RNA_GT_code =-1

        if DNA_TUMOR not in ["","."]:
            data= string.split(DNA_TUMOR,":")
            if len(string.split(data[pos],'/'))>=2:
                DNA_GT_code = int(string.split(data[pos],'/')[-1]) # the last segment
            if chrom in ["chrX","chrY"]: ##### the really stupid thing RADIA does to set GT=0 reference on chrX and Y even when there is clear evidence of altnative allele. can't believe this!
                #parse data to figure this out really stupid way to do it.
                AF_pos=-1
                data= string.split(format,":")
                for i in range(0,len(data)):
                    if data[i]== "AD":
                        AF_pos= i
                if AF_pos==-1 :
                    DNA_GT_code =-1
                elif DNA_TUMOR not in ["","."]:
                    data= string.split(DNA_TUMOR,":")
                    data = string.split(data[AF_pos],",")
                    if len(data)<2:
                        DNA_GT_code =-1
                    elif len(data)==2: # ref, alt1
                        DNA_GT_code =1
                        DNA_AD = float(data[DNA_GT_code])
                    else:
                        DNA_GT_code =1
                        DNA_AD = float(data[DNA_GT_code])
                        for i in range (2, len(data)):
                            if float(data[i]) > float(data[DNA_GT_code]):#ref, alt1, alt2, alt3:
                                DNA_GT_code = i
                                DNA_AD = float(data[DNA_GT_code])
                else:
                    DNA_GT_code =-1
        else:
            DNA_GT_code =-1

        if RNA_TUMOR not in ["","."]:
            data= string.split(RNA_TUMOR,":")
            if len(string.split(data[pos],'/'))>=2:
                RNA_GT_code = int(string.split(data[pos],'/')[-1]) # the last segment
            if chrom in ["chrX","chrY"]: ##### the really stupid thing RADIA does to set GT=0 reference on chrX and Y even when there is clear evidence of altnative allele. can't believe this!
                #parse data to figure this out myself!
                AF_pos=-1
                data= string.split(format,":")
                for i in range(0,len(data)):
                    if data[i]== "AD":
                        AF_pos=i
                if AF_pos==-1 :
                    RNA_GT_code =-1
                elif RNA_TUMOR not in ["","."]:
                    data= string.split(RNA_TUMOR,":")
                    data = string.split(data[AF_pos],",")
                    if len(data)<2:
                        RNA_GT_code =-1
                    elif len(data)==2: # ref, alt1
                        RNA_GT_code =1
                        RNA_AD = float(data[RNA_GT_code])
                    else:
                        RNA_GT_code =1
                        RNA_AD = float(data[RNA_GT_code])
                        for i in range (2, len(data)):
                            if float(data[i]) > float(data[RNA_GT_code]):#ref, alt1, alt2, alt3:
                                RNA_GT_code = i
                                RNA_AD = float(data[RNA_GT_code])
                else:
                    RNA_GT_code =-1
        else:
            RNA_GT_code =-1

        if DNA_GT_code in [-1,0] and RNA_GT_code > 0:
            return RNA_GT_code
        if RNA_GT_code in [-1,0] and DNA_GT_code >0:
            return DNA_GT_code
        if DNA_GT_code > 0 and RNA_GT_code > 0 and DNA_GT_code == RNA_GT_code:
            return DNA_GT_code
        if DNA_GT_code <= 0 and RNA_GT_code <= 0:
            return None
        if DNA_GT_code > 0 and RNA_GT_code > 0 and DNA_GT_code != RNA_GT_code and (chrom not in ["chrX","chrY"]):
            return None
        if DNA_GT_code > 0 and RNA_GT_code > 0 and DNA_GT_code != RNA_GT_code and (chrom in ["chrX","chrY"]):  # really stupid RADIA chrX and Y handling
            if RNA_AD > DNA_AD:
               return RNA_GT_code
            else:
               return DNA_GT_code

    def _parse_NORMAL_ALT_ID (self, ID, format, DNA_NORMAL, GT_code):
        #get the "ID" column in VCF
        pos=-1
        data= string.split(format,":")
        for i in range(0,len(data)):
            if data[i]== ID:
                pos=i
        if pos==-1 :
            return None

        if DNA_NORMAL not in ["","."]:
            data= string.split(DNA_NORMAL,":")
            try:
                DNA_ID_val = string.split(data[pos],",")[GT_code]
            except:
                DNA_ID_val="NA"
        else:
            DNA_ID_val="NA"

        return DNA_ID_val


    def _parse_TUMOR_ALT_ID (self, ID,format,DNA_TUMOR,RNA_TUMOR, GT_code):
        #get the "ID" column in VCF
        pos=-1
        data= string.split(format,":")
        for i in range(0,len(data)):
            if data[i]== ID:
                pos=i
        if pos==-1 :
            return None


        if DNA_TUMOR not in ["","."]:
            data= string.split(DNA_TUMOR,":")
            try:
                DNA_ID_val = string.split(data[pos],",")[GT_code]
            except:
                DNA_ID_val="NA"
        else:
            DNA_ID_val="NA"

        if RNA_TUMOR not in ["","."]:
            data= string.split(RNA_TUMOR,":")
            try:
                RNA_ID_val = string.split(data[pos],",")[GT_code]
            except:
                RNA_ID_val="NA"
        else:
            RNA_ID_val="NA"
        return [DNA_ID_val,RNA_ID_val]

    def _parse_NORMAL_SINGLE_ID (self, ID,format,DNA_NORMAL):
        #get the "ID" column in VCF
        pos=-1
        data= string.split(format,":")
        for i in range(0,len(data)):
            if data[i]== ID:
                pos=i
        if pos==-1 :
            return None

        if DNA_NORMAL not in ["","."]:
            data= string.split(DNA_NORMAL,":")
            DNA_ID_val = data[pos]
        else:
            DNA_ID_val="NA"

        return DNA_ID_val

    def _parse_TUMOR_SINGLE_ID (self, ID,format,DNA_TUMOR,RNA_TUMOR):
        #get the "ID" column in VCF
        pos=-1
        data= string.split(format,":")
        for i in range(0,len(data)):
            if data[i]== ID:
                pos=i
        if pos==-1 :
            return None

        if DNA_TUMOR not in ["","."]:
            data= string.split(DNA_TUMOR,":")
            DNA_ID_val = data[pos]
        else:
            DNA_ID_val="NA"

        if RNA_TUMOR not in ["","."]:
            data= string.split(RNA_TUMOR,":")
            RNA_ID_val = data[pos]
        else:
            RNA_ID_val="NA"
        return [DNA_ID_val,RNA_ID_val]

    def _parseEffectsPerGene(self, effectString, columnLabels, GT_code):
        if effectString =="":
            return {}
        effectPerGene = dict()
        effects = re.sub("EFF=", "", effectString).split(",")
        for thisEffect in effects:
            effectType = thisEffect.split("(")[0]
            #
            # Given a string such as
            # downstream_gene_variant(MODIFIER||3956||459|CA9|||NM_001216.2||1)
            # extract the stuff between the parens, divide it by |, and store
            # it in a dictionary indexed by the column labels given as input
            #
            effectTokens = re.sub(".+\(", "",
                                  re.sub("\)", "", thisEffect)).split("|")
            effect = dict()
            for ii in range(0,len(effectTokens)):
                effect[columnLabels[ii]] = effectTokens[ii]

            #match GT_code
            if GT_code and GT_code != int(effect["Genotype_Number"]):
                continue

            effect["effect"] = effectType
            #
            # Parse through the list of effects.  Extract the gene.
            # Save one effect per gene, choosing an arbitrary effect
            # from the most severe effect class.
            #
            thisGene = effect["Gene_Name"]
            if not effectPerGene.has_key(thisGene):
                effectPerGene[thisGene] = effect
            else:
                impactThisEffect = effect["Effect_Impact"]
                worstImpactYet = effectPerGene[thisGene]["Effect_Impact"]
                if impact[impactThisEffect] > impact[worstImpactYet]:
                    effectPerGene[thisGene] = effect
                elif impact[impactThisEffect] == impact[worstImpactYet]:
                    if effect["Amino_Acid_length"] > effectPerGene[thisGene]["Amino_Acid_length"]:
                        effectPerGene[thisGene] = effect
        return(effectPerGene)



class vcf(object):
    """This object contains the set of rows from a VCF file"""
    def __init__(self, stream):
        self._effectColumn = None
        self._rows = list()
        self._indexNextRow = 0
        for row in stream:
            if re.search("^##INFO=<ID=EFF", row):
                """Given a line of format

                ##INFO=<ID=EFF,Number=.,Type=String,Description="Pred<icted
                effects for this variant.Format: 'Effect (
                Effect_Impact | Functional_Class | Codon_Change |
                Amino_Acid_Change| Amino_Acid_length | Gene_Name |
                Transcript_BioType | Gene_Coding | Transcript_ID |
                Exon_Rank | Genotype_Number [ | ERRORS | WARNINGS ] )'">

                Parse out the ordered list of tokens as they will appear:
                ['Functional_Class', 'Codon_Change', 'Amino_Acid_Change',
                'Amino_Acid_length', 'Gene_Name', 'Transcript_BioType',
                'Gene_Coding', 'Transcript_ID', 'Exon_Rank', 'Genotype_Number'
                'ERRORS']"""
                row = re.sub("[\[\] ]", "", row)
                row = re.sub("ERRORS|WARNINGS", "ERRORS", row)
                self._effectColumn = re.split("[(|)]", row)[1:-1]
            elif not re.search("^#", row):
                # This is a row of VCF data
                newVcfRow = vcfRow(row, self._effectColumn)
                self._rows.append(newVcfRow)

    def read(self):
        return self._rows

import math

def round_sigfigs(num, sig_figs):
    """Round to specified number of sigfigs.

    >>> round_sigfigs(0, sig_figs=4)
    0
    >>> int(round_sigfigs(12345, sig_figs=2))
    12000
    >>> int(round_sigfigs(-12345, sig_figs=2))
    -12000
    >>> int(round_sigfigs(1, sig_figs=2))
    1
    >>> '{0:.3}'.format(round_sigfigs(3.1415, sig_figs=2))
    '3.1'
    >>> '{0:.3}'.format(round_sigfigs(-3.1415, sig_figs=2))
    '-3.1'
    >>> '{0:.5}'.format(round_sigfigs(0.00098765, sig_figs=2))
    '0.00099'
    >>> '{0:.6}'.format(round_sigfigs(0.00098765, sig_figs=3))
    '0.000988'
    """
    if num != 0:
        return round(num, -int(math.floor(math.log10(abs(num))) - (sig_figs - 1)))
    else:
        return 0  # Can't take the log of 0


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("ID", type=str, help="Entry for the ID column")
    parser.add_argument("output", type=str, help="outputfile")
    args = parser.parse_args()

    myVcf = vcf(sys.stdin)

    fout =open(args.ID,'w')
    total=0
    good=0
    for row in myVcf.read():
        total =total+1
        if row.alt =="NA":
            continue ######## bad calls in the VCF
        good = good+1
        if str(row.DNA_AF) not in ["NA",""]:
            row.DNA_AF= round_sigfigs(float(row.DNA_AF),3)
        else:
            row.DNA_AF=""
        if str(row.RNA_AF) not in ["NA",""]:
            row.RNA_AF= round_sigfigs(float(row.RNA_AF),3)
        else:
            row.RNA_AF=""
        if str(row.NORMAL_AF) not in ["NA",""]:
            row.NORMAL_AF= round_sigfigs(float(row.NORMAL_AF),3)
        else:
            row.NORMAL_AF=""
        if len(row.effectPerGene)!=0:
            for gene in row.effectPerGene.keys():
                AA_Change = row.effectPerGene[gene]["Amino_Acid_Change"]
                if AA_Change !="" and AA_Change[:2]!="p.":
                    AA_Change="p."+AA_Change
                fout.write(string.join([args.ID, row.chr, str(row.start),
                                        str(row.end), row.reference, row.alt,
                                        gene,row.effectPerGene[gene]["effect"],
                                        str(row.DNA_AF), str(row.RNA_AF),AA_Change,
                                        str(row.NORMAL_AF)],"\t")+"\n")
        else:
            gene =""
            AA_Change=""
            effect =""
            fout.write(string.join([args.ID, row.chr, str(row.start),
                                    str(row.end), row.reference, row.alt,
                                    gene,effect,
                                    str(row.DNA_AF), str(row.RNA_AF),AA_Change,
                                    str(row.NORMAL_AF)],"\t")+"\n")
    fout.close()

    if float(good)/float(total)>0.9:
        os.system("cat "+args.ID+" >> "+args.output)
        os.system("rm -f "+args.ID)
    else:
        os.system("rm -f "+args.ID)
        print "throw out due to "+ str(1- float(good)/float(total))+ " calls in vcf dose not make sense", args.ID

if __name__ == '__main__':
        main()


