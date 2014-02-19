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
        self.chr = tokens[0]
        self.start = int(tokens[1])
        self.strand = tokens[2]
        self.reference = tokens[3]
        self.end = self.start + len(self.reference) - 1

        # this could be all specific to RADIA output
        format =tokens[8] 
        DNA_TUMOR =tokens[10] 
        RNA_TUMOR = tokens[11]
        
        #get the ALT base code (in VCF: GT=0,1,2?, 1 and 2 are acceptable)
        GT_code = self._findGTCode(format, DNA_TUMOR, RNA_TUMOR)
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
            val =self._parseDNA_TUMOR_ALT_ID(ID,format,DNA_TUMOR, RNA_TUMOR)
            if val != None:
                self.DNA_AF, self.RNA_AF= val
            else: # returned None
                #print "AF error", row
                self.DNA_AF, self.RNA_AF= "NA","NA"

        for thisSubToken in tokens[7].split(";"):
            if re.search("^EFF=", thisSubToken):
                effectsString = thisSubToken
        self.effectPerGene = self._parseEffectsPerGene(effectsString, columnLabels)

    def _findGTCode(self, format, DNA_TUMOR, RNA_TUMOR):
        pos=-1 
        data= string.split(format,":")
        for i in range(0,len(data)):
            if data[i]== "GT":
                pos=i
                break
        if pos ==-1:
            return None

        if DNA_TUMOR not in ["","."]:
            data= string.split(DNA_TUMOR,":")
            DNA_GT_code = int(string.split(data[pos],'/')[-1]) # the last segment
        else:
            DNA_GT_code =-1

        if RNA_TUMOR not in ["","."]:
            data= string.split(RNA_TUMOR,":")
            RNA_GT_code = int(string.split(data[pos],'/')[-1]) # the last segment
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
        if DNA_GT_code > 0 and RNA_GT_code > 0 and DNA_GT_code != RNA_GT_code:
            return None
        
    def _parseDNA_TUMOR_ALT_ID (self, ID,format,DNA_TUMOR,RNA_TUMOR):
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
            DNA_ID_val = string.split(data[pos],",")[-1]  # the last segment is for the ALT allele
        else:
            DNA_ID_val="NA"
            
        if RNA_TUMOR not in ["","."]:
            data= string.split(RNA_TUMOR,":")
            RNA_ID_val = string.split(data[pos],",")[-1]  # the last segment is for the ALT allele
        else:
            RNA_ID_val ="NA"
        return [DNA_ID_val,RNA_ID_val]
    
    def _parseEffectsPerGene(self, effectString, columnLabels):
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
    args = parser.parse_args()

    myVcf = vcf(sys.stdin)
    for row in myVcf.read():
        for gene in row.effectPerGene.keys():
            print "%s\t%s\t%d\t%d\t%s\t%s\t%s\t%s\t%s\t%s" % (args.ID, row.chr, row.start,
                                                  row.end, gene, row.reference, row.alt, 
                                                  row.effectPerGene[gene]["effect"], row.DNA_AF, row.RNA_AF)
            
if __name__ == '__main__':
        main()

        
    
    
                
    
