import string, os, sys

def findAve(DNA_AF):
    n=0
    total =0.0
    for item in DNA_AF:
        if item in ["","NA","N/A","NULL", "null","Null"]:
            continue
        n= n+1
        total = total +  float(item)
    if n==0:
        return "NA"
    else:
        return str(total/n)


if len (sys.argv[:])!=3:
    print "python removeDuplicate.py inputXena outputXena"
    sys.exit()

fin = open(sys.argv[1],'r')
fout = open(sys.argv[2],'w')


dic ={}
for line in fin.readlines():
    line = string.strip(line)
    if line[:6] =="sample":
        fout.write(line+"\n")
        continue
    data = string.split(line,'\t')
    sample = data[0]
    chrom = data[1]
    start = data[2]
    end = data[3]
    gene = data[4]
    ID = sample+chrom+start+end+gene
    if dic.has_key(ID):
        dic[ID].append(data)
    else:
        dic[ID]=[data]
fin.close()

for ID in dic.keys():
    if len(dic[ID])==1:
        fout.write(string.join(dic[ID][0],"\t")+"\n")
    else:
        effects=[]
        alts=[]
        for i in range (0,len(dic[ID])):
            e=  dic[ID][i][7]
            alt = dic[ID][i][6]
            if e not in effects:
                effects.append(e)
            if alt not in alts:
                alts.append(alt)
        if len(effects)!=1 or len(alts)!=1 :
            print dic[ID],effects, alts
            continue ###################### 
        
        #compute average DAN_AF and RNA_AF
        DNA_AF=[]
        RNA_AF=[]
        for i in range (0,len(dic[ID])):
            d = dic[ID][i][8]
            r = dic[ID][i][9]
            DNA_AF.append(d)
            RNA_AF.append(r)
        
        ave_DNA_AF = findAve(DNA_AF)
        ave_RNA_AF = findAve(RNA_AF)
        dic[ID][0][8]= ave_DNA_AF
        dic[ID][0][9]= ave_RNA_AF
        fout.write(string.join(dic[ID][0],"\t")+"\n")
fout.close()

