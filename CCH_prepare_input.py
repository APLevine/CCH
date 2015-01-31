####
##  Convert TPED/TFAM 12 (Plink) format to SNP format
##  A. P. Levine
##  March 2014
####

## Execution syntax:
## python CCH_prepare_input.py input.tped input.tfam > output.snp

import sys

ped=open(sys.argv[1],"r") #tped file
fam=open(sys.argv[2],"r") #tfam file

id=[]
for line in fam:
    split=line.rstrip("\n").split()
    id.append(split[1])

print "Chr\tPos_(bp)\tPosition\tSNP_ID\tAllele_A\tAllele_B\t"+"\t".join(id)

g=""
g2=""
for line in ped:
    split=line.rstrip("\n").split()
    count = 0
    for i in split[4:]:
        count = count + 1
        if i=="1" or i=="A":
            j="A"
        if i=="2" or i=="B":
            j="B"
        if i=="0":
            j="0"
        if count%2==0:
            g=g+j
        if count%2!=0:
            g=g+"\t"+j
    for i in g.split("\t"):
        if i=="AB":
            g2=g2+"\t"+"AB"
        if i=="BA":
            g2=g2+"\t"+"AB"
        if i=="AA":
            g2=g2+"\t"+"AA"
        if i=="BB":
            g2=g2+"\t"+"BB"
        if i=="00":
            g2=g2+"\t"+"00"
    print split[0]+"\t"+split[3]+"\t"+split[2]+"\t"+split[1]+"\t?\t?"+ g2
    g=""
    g2=""

