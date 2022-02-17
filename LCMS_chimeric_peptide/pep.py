#! /usr/bin/env python3

'''
Author: Ajeet Singh
Email: singh.ajeet@nipgr.ac.in

This script is designed to get chimeric junction amino acid containing trypisn digested peptides
from 3'Frame translated chimeric RNA with 200bp from each participating gene from their break-position.
These chimeric RNAs were detected from paired-end RNA-Seq using tools like Star-Fusion, Trinity-Fusion,
Squid, Mapsplice, FusionBloom.

Usage:

    python3 pep.py SAMPLE

"SAMPLE" is directory, should be present in sample working directory containing all final results files
from above mentioned tools. Provide "SAMPLE" only, without "/", to "pep.py". It gives output "SAMPLE_final_pep.fa"
'''

import sys
import bisect
import os
import numpy as np
import subprocess
from Bio import SeqIO

srr_id=sys.argv[1]

starf        = srr_id+"/star-fusion.fusion_predictions.abridged.coding_effect.tsv"
trinityf     = srr_id+"/TrinityFusion-UC.fusion_predictions.tsv"
ericf        = srr_id+"/fusion_eric_arabidopsis.results.filtered.tsv"
mapsplicef   = srr_id+"/fusions_candidates.txt"
squidf       = srr_id+"/squid_out_sv.txt"
fusionbloomf = srr_id+"/sv.bedpe"

removf=['rm','-f', srr_id+"_final_pep.fa"]
try:
    subprocess.run(removf, check=True)
except subprocess.CalledProcessError as e:
    print("\nError!\n")
    sys.exit(e.stderr)

posi=[]
def charposition(mstring , mchar):
    for n in range(len(mstring)):
        if mstring[n] == mchar:
            posi.append(n)
    return posi


def bedprint(g1a,g1b,g1c,g1d,g2a,g2b,g2c,g2d):
    if g1d=="+":
        bedstr1=g1b, int(g1c)-200, g1c, g1a+"_"+g2a+"::"+srr_id, 0, g1d
    elif g1d=="-":
        bedstr1=g1b, g1c, int(g1c)+200, g1a+"_"+g2a+"::"+srr_id, 0, g1d

    if g2d=="+":
        bedstr2=g2b, g2c, int(g2c)+200, g1a+"_"+g2a+"::"+srr_id, 0, g2d
    elif g2d=="-":
        bedstr2=g2b, int(g2c)-200, g2c, g1a+"_"+g2a+"::"+srr_id, 0, g2d

    if g1d==".":
        bedstr1=g1b, int(g1c)-200, g1c, g1a+"_"+g2a+"::"+srr_id, 0, g1d
    if g2d==".":
        bedstr2=g2b, g2c, int(g2c)+200, g1a+"_"+g2a+"::"+srr_id, 0, g2d

    tuple_list=[bedstr1,bedstr2]
    with open(srr_id+'.bed', 'w') as f:
        for tuple in tuple_list:
            f.write('%s\t%s\t%s\t%s\t%s\t%s\n' % tuple)
    
    genome="/home/ajeet/Projects/PFusion/LCMS_chimeric_peptide/genome.fa"
    barg = ['bedtools', 'getfasta', '-fi', genome, '-bed', srr_id+'.bed', '-s', '-name']
    try:
        fa = open(srr_id+".fa", "w")
        subprocess.run(barg, check=True, stdout=fa)
        fa.close()
    except subprocess.CalledProcessError as e:
        print("\nError!\n")
        sys.exit(e.stderr)

    #concate the both sequences side by side
    heade=[]
    seq=[]
    with open(srr_id+".fa", 'r') as f:
        lines=f.readlines()
        for x in lines:
            if ">" in x:
                heade.extend(x.lstrip(">").rstrip("\n").split())
            if not ">" in x:
                seq.extend(x.rstrip("\n").split())
        mheade = ''.join(heade)
        mseq = ''.join(seq)
        with open(srr_id+".concat.fa", 'w') as f:
            f.write(">%s\n" % mheade)
            f.write("%s\n" % mseq)

    # Translate the fasta seq 
    transeq="/home/ajeet/Software/EMBOSS-6.6.0/bin/transeq"
    targ = [transeq, '-sequence', srr_id+".concat.fa", '-outseq',  srr_id+"_pep.fa", '-frame', 'F', '-trim', 'Y']
    try:
        er=open("error",'w')
        subprocess.run(targ, check=True, stderr=er)
        er.close()
    except subprocess.CalledProcessError as e:
        print("\nError!\n")
        sys.exit(e.stderr)

    #convert fasta (multi-line) to single-line tabular fasta and get sequence region with aa from junction codon
    mystr=""
    for record in SeqIO.parse(srr_id+"_pep.fa", 'fasta'):
        list='>{}\t{}'.format(record.description, record.seq)
        if not "*" in list.split('\t')[1][64:67]:
            mystr=list.split('\t')[1][64:67]
        else:
            mystr="XX"

        for seq in list.split('\t')[1].split():
            if "*" in seq:
                posi = charposition(seq, '*')
                i = bisect.bisect_left(posi , 65)
                mend = posi[i]
                mstart = posi[i - 1]
                chunke = seq[mstart +1 :mend]
                with open(srr_id+"_sel_pep.fa", 'w') as f:
                    f.write("%s\n" % heade)
                    f.write("%s\n" % chunke)

            if not "*" in seq:
                if seq.endswith(mystr, 64, 67) == True:
                    with open(srr_id+"_sel_pep.fa", 'w') as f:
                        f.write("%s\n" % heade)
                        f.write("%s\n" % seq)

    darg = ['/home/ajeet/Software/EMBOSS-6.6.0/bin/pepdigest', '-seqall', srr_id+"_pep.fa", '-outfile', srr_id+".pepdigest", '-menu', '1', '-mono', 'no']
    try:
        er=open("error",'w')
        subprocess.run(darg, check=True, stderr=er)
        er.close()
    except subprocess.CalledProcessError as e:
        print("\nError!\n")
        sys.exit(e.stderr)
                    
    #fetch trypsin digested peptides with chimeric junction amino acids
    with open(srr_id+".pepdigest", 'r') as f:
        for i in f.readlines():
            if not i.strip():
                continue
            if i:
                if not i.strip().startswith(("#", "Start")):
                    if len(i.strip().split(' ')[-1]) >= 6:
                        for seq in i.strip().split(' ')[-1].split():
                            if mystr in seq:
                                with open(srr_id+"_final_pep.fa", 'a') as f:
                                    f.write(">%s\n" % mheade)
                                    f.write("%s\n" % seq)
                                
    
gene1=[]
chr_gene1=[]
bp_gene1=[]
strand_gene1=[]
gene2=[]
chr_gene2=[]
bp_gene2=[]
strand_gene2=[]

def starfBED(starf):
    if (os.path.isfile(starf)):
        #print ("STAR-Fusion results exits for "+srr_id)

        with open(starf,'r') as f:
            next(f)
            lines=f.readlines()
            for x in lines:
                gene1=x.split('\t')[6]
                chr_gene1=x.split("\t")[7].split(":")[0]
                bp_gene1=x.split("\t")[7].split(":")[1]
                strand_gene1=x.split("\t")[7].split(":")[2]

                gene2=x.split('\t')[8]
                chr_gene2=x.split("\t")[9].split(":")[0]
                bp_gene2=x.split("\t")[9].split(":")[1]
                strand_gene2=x.split("\t")[9].split(":")[2]

                bedprint(gene1,chr_gene1,bp_gene1,strand_gene1,gene2,chr_gene2,bp_gene2,strand_gene2)

def trinityfBED(trinityf):
    if (os.path.isfile(trinityf)):
        #print ("Trinity-Fusion results exits for "+srr_id)
        with open(trinityf,'r') as f:
            next(f)
            lines=f.readlines()
            for x in lines:
                gene1=x.split('\t')[5]
                chr_gene1=x.split("\t")[6].split(":")[0]
                bp_gene1=x.split("\t")[6].split(":")[1]
                strand_gene1="."

                gene2=x.split('\t')[7]
                chr_gene2=x.split("\t")[8].split(":")[0]
                bp_gene2=x.split("\t")[8].split(":")[1]
                strand_gene2="."

                bedprint(gene1,chr_gene1,bp_gene1,strand_gene1,gene2,chr_gene2,bp_gene2,strand_gene2)

def ericfBED(ericf):
    if (os.path.isfile(ericf)):
        #print ("EricScripts-Plants results exits for "+srr_id)
        with open(ericf,'r') as f:
            next(f)
            lines=f.readlines()
            for x in lines:
                gene1=x.split('\t')[8]
                chr_gene1=x.split("\t")[2]
                bp_gene1=x.split("\t")[3]
                strand_gene1=x.split("\t")[4]

                gene2=x.split('\t')[9]
                chr_gene2=x.split("\t")[5]
                bp_gene2=x.split("\t")[6]
                strand_gene2=x.split("\t")[7]

                bedprint(gene1,chr_gene1,bp_gene1,strand_gene1,gene2,chr_gene2,bp_gene2,strand_gene2)

def mapsplicefBED(mapsplicef):
    if (os.path.isfile(mapsplicef)):
        #print ("MapSplice2 results exits for "+srr_id)
        with open(mapsplicef,'r') as f:
            next(f)
            lines=f.readlines()
            for x in lines:
                gene1='gene1'
                chr_gene1=x.split("\t")[0].split("~")[0]
                bp_gene1=x.split("\t")[1]
                strand_gene1=x.split("\t")[5][0]

                gene2='gene2'
                chr_gene2=x.split("\t")[0].split("~")[1]
                bp_gene2=x.split("\t")[2]
                strand_gene2=x.split("\t")[5][1]

                bedprint(gene1,chr_gene1,bp_gene1,strand_gene1,gene2,chr_gene2,bp_gene2,strand_gene2)

def squidfBED(squidf):
    if (os.path.isfile(squidf)):
        #print ("SQIUD results exits for "+srr_id)
        with open(squidf,'r') as f:
            next(f)
            lines=f.readlines()
            for x in lines:
                gene1='gene1'
                chr_gene1=x.split("\t")[0]
                chk_str1=x.split("\t")[8]
                if chk_str1=="-":
                    bp_gene1=x.split("\t")[1]
                elif chk_str1=="+":
                    bp_gene1=x.split("\t")[2]
                strand_gene1=x.split("\t")[8]

                gene2='gene2'
                chr_gene2=x.split("\t")[3]
                chk_str2=x.split("\t")[9]
                if chk_str2=="-":
                    bp_gene2=x.split("\t")[4]
                elif chk_str2=="+":
                    bp_gene2=x.split("\t")[5]
                strand_gene2=x.split("\t")[9]

                bedprint(gene1,chr_gene1,bp_gene1,strand_gene1,gene2,chr_gene2,bp_gene2,strand_gene2)

def fusionbloomfBED(fusionbloomf):
    if (os.path.isfile(fusionbloomf)):
        #print ("FusionBloom results exits for "+srr_id)
        with open(fusionbloomf,'r') as f:
            lines=f.readlines()[4:]
            for x in lines:
                gene1=x.split('\t')[14]
                chr_gene1=x.split("\t")[0]
                bp_gene1=x.split("\t")[2]
                strand_gene1=x.split("\t")[8]

                gene2=x.split('\t')[19]
                chr_gene2=x.split("\t")[3]
                bp_gene2=x.split("\t")[4]
                strand_gene2=x.split("\t")[9]

                bedprint(gene1,chr_gene1,bp_gene1,strand_gene1,gene2,chr_gene2,bp_gene2,strand_gene2)

if __name__ == '__main__':
    starfBED(starf)
    trinityfBED(trinityf)
    ericfBED(ericf)
    mapsplicefBED(mapsplicef)
    squidfBED(squidf)
    fusionbloomfBED(fusionbloomf)


removef=['rm','-f', srr_id+'.bed', srr_id+'.fa', srr_id+".concat.fa", srr_id+".pepdigest", srr_id+"_pep.fa", srr_id+"_sel_pep.fa", 'error']
try:
    subprocess.run(removef, check=True)
except subprocess.CalledProcessError as e:
    print("\nError!\n")
    sys.exit(e.stderr)

