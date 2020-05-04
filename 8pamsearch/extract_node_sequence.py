#!/bin/bash
# Tue Jul 30 10:52:02 CEST 2019
# Made by L-F-S
# At the University Of Trento, Italy

import os
import sys
import numpy as np
import pandas as pd

sys.path.insert(0, '/home/lorenzo.signorini/utils/')
import filename_discrepancies
from Bio import SeqIO

dataset="CM_madagascar" #sys.argv[1] #TODO switch test dataset
outdir="/shares/CIBIO-Storage/CM/scratch/tmp_projects/signorini_cas/Nstep" #TODO add nstepdir
datadir="/shares/CIBIO-Storage/CM/scratch/tmp_projects/signorini_cas" #TODO add nstepdir
annodir="/shares/CIBIO-Storage/CM/scratch/tmp_projects/signorini_cas/2casanno/crisprcasanno"
allcontsdir="/shares/CIBIO-Storage/CM/scratch/tmp_projects/epasolli_metagenomicassembly"
contig=sys.argv[1]#"NODE_59_length_82873_cov_9.4706"
sample="A01_02_1FE"

s3dataset=dataset #TODO filename discrepancies!
os.chdir(allcontsdir+"/"+dataset+"/"+sample) #TODO filename discrepancies
allconts=open("contigs_filtered.fasta", "r")
contseq=""
for record in SeqIO.parse("contigs_filtered.fasta", "fasta"):
    if record.id.startswith(contig):
        contseq=record
allconts.close()
#write file
SeqIO.write(contseq,datadir+"/4exploratoryanalyses/pamsearch/test/"+contig+"_putative_viral.fasta","fasta")
#f=open(datadir+"/4exploratoryanalyses/pamsearch/test/"+contig+"_putative_protospacer.fasta")  #TODO TEST folder!
#f.write(">"+contig+"\n"+contseq)
#f.close()
