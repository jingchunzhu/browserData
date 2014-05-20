#!/usr/bin/env python

import argparse
import re
import sys
import string

impact = {"MODIFIER":0, "LOW":1, "MODERATE":2, "HIGH":3}

class vcfRow(object):
    """This object contains an individual row of VCF output"""
    def __init__(self, row, columnLabels):
        tokens = row[:-1].split("\t")
        #chromosome notation
        chrom = tokens[0]
        if string.upper(chrom[0:3])== "CHR": 
            chrom[0:3]="chr"
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
        self.end = self.start + len(self.reference) - 1

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
            DNA_TUMOR =tokens[10] 
            if len(tokens)>11:
                RNA_TUMOR = tokens[11]
            else:
                RNA_TUMOR=''
            
            #get the ALT base code (in VCF: GT=0,1,2?, 1 and 2 are acceptable)
            GT_code = self._findGTCode(self.chr, format, DNA_TUMOR, RNA_TUMOR, self.start)
            if GT_code !=None:
                self.alt = string.split(tokens[4],",")[GT_code-1]
            else:
                #print "GT error", row
                self.alt = "NA"
        
            # AF allel frequency # this is all specific to RADIA output
            if GT_code == None:
                self.DNA_AF, self.RNA_AF= "NA","NA"
            else:
                ID="AF"
                val =self._parse_TUMOR_ALT_ID(ID,format,DNA_TUMOR, RNA_TUMOR, GT_code)
                if val != None:
                    self.DNA_AF, self.RNA_AF= val
                else: # returned None
                    #print "AF error", row
                    self.DNA_AF, self.RNA_AF= "NA","NA"

        self.effectPerGene = self._parseEffectsPerGene(effectsString, columnLabels, GT_code)

    def _findGTCode(self, chrom, format, DNA_TUMOR, RNA_TUMOR, start):
        pos=-1 
        data= string.split(format,":")
        for i in range(0,len(data)):
            if data[i]== "GT":
                pos=i
                break
        if pos ==-1:
            return None

        if chrom in ["chrX","chrY"]:
            DNA_AD=-1
            RNA_AD=-1

        if DNA_TUMOR not in ["","."]:
            data= string.split(DNA_TUMOR,":")
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
                    if len(data)==2: # ref, alt1
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
                    if len(data)==2: # ref, alt1
                        RNA_GT_code =1
                        RNA_AD = float(data[RNA_GT_code])
                    else:
                        RNA_GT_code =1
                        RNA_AD = float(data[RNA_GT_code])
                        print data
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
            DNA_ID_val = string.split(data[pos],",")[GT_code]  
        else:
            DNA_ID_val="NA"
            
        if RNA_TUMOR not in ["","."]:
            data= string.split(RNA_TUMOR,":")
            RNA_ID_val = string.split(data[pos],",")[GT_code]  
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


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("ID", type=str, help="Entry for the ID column")
    parser.add_argument("output", type=str, help="outputfile")
    args = parser.parse_args()

    myVcf = vcf(sys.stdin)
    fout =open(args.output,'a')
    for row in myVcf.read():
        if len(row.effectPerGene)!=0:
            for gene in row.effectPerGene.keys():
                AA_Change = row.effectPerGene[gene]["Amino_Acid_Change"]
                fout.write(string.join([args.ID, row.chr, str(row.start),
                                        str(row.end), gene, row.reference, row.alt, 
                                        row.effectPerGene[gene]["effect"], str(row.DNA_AF), str(row.RNA_AF),AA_Change],"\t")+"\n")
        else:
            fout.write(string.join([args.ID, row.chr, str(row.start),
                                    str(row.end), "", row.reference, row.alt, 
                                    "", str(row.DNA_AF), str(row.RNA_AF),""],"\t")+"\n")
    fout.close()
if __name__ == '__main__':
        main()

        
    
    
                
    
