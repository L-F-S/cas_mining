#!/usr/bin/env python
# coding: utf-8

# # find tracerRNA
# 

# In[1]:


# First created 08-10-2019
# Made by L-F-S
# At the University Of Trento, Italy
#
# USAGE:
    # python find_tracrRNA_of_id <CAS9_ID>

import os
import sys
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Alphabet import generic_dna, generic_protein, RNAAlphabet
from Bio import pairwise2


sys.path.insert(0, '/home/lorenzo.signorini/cas_mining/utils/')
import filename_discrepancies

feature="Cas9" #WARNING!!! CHANGE THIS!!
#outpath="/shares/CIBIO-Storage/CM/scratch/tmp_projects/signorini_cas/5caslocitable/out/"+feature+"/"
datadir="/shares/CIBIO-Storage/CM/scratch/tmp_projects/signorini_cas/"



# In[2]:


# Import data
cas_dataset=pd.read_csv("/shares/CIBIO-Storage/CM/scratch/tmp_projects/signorini_cas/5caslocitable/known_"+feature+"_variants_table.csv", index_col=0)
cas_dataset.columns


# In[3]:


tabellazza=pd.read_csv(datadir+"/3tabellazza/crisprcas_hits_table.csv", index_col=0)
tabellazza.columns


# In[4]:


seqid =sys.argv[1] #"CJFLGLGO_01367"
print("SEQ ID: ", seqid)
print(cas_dataset[cas_dataset["Seq ID"]==seqid]["Contig"].iloc[0])
print(cas_dataset[cas_dataset["Seq ID"]==seqid]["Study"].iloc[0])
print(cas_dataset[cas_dataset["Seq ID"]==seqid]["Genome Name"].iloc[0])
print(len(cas_dataset[cas_dataset["Seq ID"]==seqid]["Seq"].iloc[0]))
print(cas_dataset[cas_dataset["Seq ID"]==seqid]["minced_CRISPR"].iloc[0])
print("----------------------------------------------------------------------")


# In[5]:


spacers=[]
repeats=[]
repeat_start_pos=[]


cas9_aa=cas_dataset[cas_dataset["Seq ID"]==seqid]["Seq"].iloc[0]
dataset=cas_dataset[cas_dataset["Seq ID"]==seqid]["Study"].iloc[0]
genomename=cas_dataset[cas_dataset["Seq ID"]==seqid]["Genome Name"].iloc[0]
SGB=cas_dataset[cas_dataset["Seq ID"]==seqid]["SGB ID"].iloc[0]
contigname=cas_dataset[cas_dataset["Seq ID"]==seqid]["Contig"].iloc[0]
print("\n-----------------------------------------------------\n"+feature+" id:\t"+seqid+"\nSGB:\t"+str(SGB)+"\nGenome Name:\t",      genomename+"\nContig:\t"+cas_dataset[cas_dataset["Seq ID"]==seqid]["Contig"].iloc[0],"\n-----------------------------------------------------\n")

######################################################################TODO THIS MUST BE INSIDE filename_discrepancies
# take into account filename discrepancies
if dataset.startswith("ZeeviD"):  #TODO guarda se funziona quando finisci di girare minced_runner
         s3_dataset="ZeeviD_2015"
         #genomename=filename_discrepancies.change_to_megahit(genomename)
         s3_genomename=genomename
         try:
             genomename=genomename.replace("ZeeviD_2015", "ZeeviD_2015_B")
             prokkafile_ofsample=open(annodir+"/justminced/ZeeviD_2015_B/"+genomename+".crisprcas.gff.minced")
             prokkafile_ofsample.close()
             dataset="ZeeviD_2015_B"
         except:
             genomename=genomename.replace("ZeeviD_2015_B", "ZeeviD_2015_A")
             dataset="ZeeviD_2015_A"
else:
    s3_dataset=dataset
    s3_genomename=genomename
    dataset=filename_discrepancies.s3(s3_dataset,r=True) # get working dataset name from s3 dataset name
    genomename=genomename.replace(s3_dataset, dataset)

if filename_discrepancies.dataset_has_megahit(dataset,genomename):
    genomename=filename_discrepancies.change_to_megahit(genomename)

samplename = genomename.split("__")[1]
########################################################################

# get CRISPR spacer and repeat sequence

mincedCRISPRfilename=datadir+"1crisprsearch/out/"+s3_dataset+"/"+s3_genomename+".fa.minced.out"
f=open(mincedCRISPRfilename, "r")
for line in f.readlines():
    print(line)
    if line[0]== "1" or line[0]== "2" or line[0]== "3" or line[0]== "4" or line[0]== "5" or     line[0]== "6" or line[0]== "7" or line[0]== "8" or line[0]== "9" or line[0]== "0":

        repeats.append(line.split("\t")[2])

        repeat_start_pos.append(line.split("\t")[0])

        if not line.split("\t")[3]=="\n":
            spacers.append(line.split("\t")[3])
f.close()

# get whole locus sequence todo this could become a function of its own, dentro utils, gli dai genomename della cas e/o seqid della cas (anzi famo solo seqid, che magari ci sn doppioni in genomi) e ti dà tt il locus

# print start and stop of CRISPR
unostranosplit=cas_dataset[cas_dataset["Seq ID"]==seqid]["minced_CRISPR"].iloc[0].split("\', \'")
crispr_start, crispr_stop = int(unostranosplit[1]), int(unostranosplit[2].strip("\')]\""))

# > get start and stop position of all Cass. questo devi metterlo nella tabella di cas non ci sono santi.

prokka_anno=tabellazza[tabellazza["Genome Name"]==s3_genomename]["prokka_cas"].iloc[0]
# so già che all casses inside here are on the same contig which is 'contig', xke le ho estratte per essere così

prokkaslist=prokka_anno.split(">")
prokkaslist.pop(0)
# get contigname of current effectorcas   
cas_position={}
for casanno in prokkaslist:
    unostranosplit=casanno.split("Name=")
    casnumber=unostranosplit[1][3] #because we know that we are loooking at 'good' loci, this is just going to be 9, 1 or 2
    unaltrostranosplit=casanno.split("\t")
    cas_start=unaltrostranosplit[3]
    cas_stop=unaltrostranosplit[4]
    cas_position["Cas"+casnumber]=[int(cas_start),int(cas_stop)]

genomefilename= "/shares/CIBIO-Storage/CM/scratch/tmp_projects/epasolli_darkmatter/allcontigs/ALLreconstructedgenomes/"+str(SGB)+"/"+genomename+".fa"
print("Retrieving locus")


# access contig:

for record in SeqIO.parse(genomefilename,"fasta"):  #c'è modo di non fare questo 'ciclo'?
    if record.id.startswith(contigname):
        #test se effettivamente la seq di cas è dove voglio io

        SeqIO.write(record, "/shares/CIBIO-Storage/CM/scratch/tmp_projects/signorini_cas/7tracrRNA/tempblastdb", "fasta")
        print("CRISPR si trova", crispr_start, crispr_stop)
        print("Cas9 si trova", cas_position["Cas9"])
        print("Cas2 si trova", cas_position["Cas2"])
        print("Cas1 si trova", cas_position["Cas1"])

        #estrai la posizione di cas9 dal contig, e traducila
        cas9start=[int(x) for x in cas_position["Cas9"]][0]
        cas9stop=[int(x) for x in cas_position["Cas9"]][1]
        cas9_nnseq=record.seq[cas9start-1:cas9stop-3] #gff should have a 1-based positional annotation (ma sto andando a occhio finchè non sono uguali, è già un oggetto Bio.Seq

        my_translated_cas9=cas9_nnseq.transcribe().translate()

  
        print(len(cas9_nnseq)/3,len(cas9_aa))


        if my_translated_cas9==cas9_aa:
            print("Cas locus on plus strand")
            plus=True
        if not my_translated_cas9==cas9_aa:
            #proviamo col reverse complement
            cas9_nnseq=record.seq[cas9start+2:cas9stop]
            my_translated_cas9=cas9_nnseq.reverse_complement().transcribe().translate()
            print("Cas locus on minus strand")
        print("Is Cas9 translated from the genome equal to the orignal annotation?", my_translated_cas9==cas9_aa)
        print(my_translated_cas9)

        #extact a locus, 200 bp a monte e a valle dei cosi che ho:
        starts, stops=[startstop[0]for startstop in cas_position.values()], [startstop[1]for startstop in cas_position.values() ]
        minimo=min([crispr_start]+starts)
        print("start: ", minimo)
        massimo=max([crispr_stop]+stops)
        print("end: ", massimo)

        #opzione 1 blast
      #1. salva np.uniq(repeats) in un fasta file: build temporary query file for blastn search

        print("let's now blast the repeat against  the genome")
        blast_folder="/shares/CIBIO-Storage/CM/scratch/tmp_projects/signorini_cas/7tracrRNA/"
        os.chdir(blast_folder)
        tempfile=open("temp_repeat_seq", "w")
        tempfile.close()
        tempfile=open("temp_repeat_seq", "a")
        for n, repeat in enumerate([Seq(sequence) for sequence in np.unique(repeats)]):
            tempfile.write(">rpt"+str(n+1)+"|"+contigname+"|"+genomename+"|"+seqid+"\n"+str(repeat)+"\n") 
        tempfile.close()

        #2. Blast query file against db
        print("Running blastn of query file agianst all contigs")

        #3. make db file
        dbfile="/shares/CIBIO-Storage/CM/scratch/tmp_projects/signorini_cas/7tracrRNA/tempblastdb"
        os.system("makeblastdb -in "+dbfile+" -parse_seqids  -dbtype nucl")



        blastoutfile="/shares/CIBIO-Storage/CM/scratch/tmp_projects/signorini_cas/7tracrRNA/"        +str(SGB)+"__"+seqid+"__"+genomename+"__"+feature+".trcr.blastout"
        blastn_command = "blastn -out "+blastoutfile+" -outfmt \"6  qseqid sseqid pident qlen length mismatch gapopen qseq sseq sstart send evalue sstrand\" -query temp_repeat_seq -db "+dbfile+" -evalue 0.001 -word_size 11"
        print(blastn_command)
        os.system(blastn_command)
       # os.system("rm /shares/CIBIO-Storage/CM/scratch/tmp_projects/signorini_cas/7tracerRNA/temp*")


# In[6]:


#4. parse blast output: only print output if it is NOT one of the old repeats.
print("\n-----------------------------------------------------\nParsing BLAST output, looking for hit which are NOT the same repeats\ni.e. not starting with:", repeat_start_pos,":\n-----------------------------------------------------\n")
f=open(blastoutfile,"r")
for line in f.readlines():
    start_pos=line.split("\t")[9]
    if start_pos not in repeat_start_pos:
        if str(int(start_pos)+1) not in repeat_start_pos:
            if str(int(start_pos)-1) not in repeat_start_pos:
                print(line) 


# ## extract sequences

# In[7]:

input_matched_rep_start=int(input"input matched repetition start position (i.e. actual lowest number)"))
input_matched_rep_stop=int(input"input matched repetition stop position (i.e. actual lowest number)"))
matched_rep_start=input_matched_rep_start
matched_rep_stop=input_matched_rep_stop


# In[26]:


# this is the sequence matching the repeat:
    matched_rep=Seq(str(input("input matched repeat sequence (from lines above): "))#Seq("TTGTAGTAACCCAATTGTTCTTGGTATGTTAT")

#rpts which match:
rpt=Seq(str(input("Input original repeat sequence (chose your favourite from 'minced' file lines above): ")))
rpt1=Seq("ATTATAGTTTCTTAATTTTTCTTGGTATGTTATAAT")
rpt2=Seq("ATTGTAGTTCCCTAATTTTTCTTGGTATATTATAAT")
rpt3=Seq("ATTGTAGTTCCCTAATTTTTCTTGGTATGTTATAAT")

# since repeat is on the + strand, and this is all on the + strand, tracrRNA should be its reverse complement:
antirepeat=matched_rep.reverse_complement()
#this should go frombo to bo  = x nt  (x+1)
print(antirepeat)
len(antirepeat)


# In[13]:


# print terminators

f=open(blast_folder+str(SGB)+"__"+seqid+"__"+genomename+".terminators")
for line in f.readlines():
    if line.strip().startswith("4"):
        print(line)
f.close()


# In[15]:


# this is the rho-independent trnascription terminator, going from 14719 or 14722 to bo (+1 perché conti il primo tipo), 
input_term_end=int(input("Input terminator end position (start position  of the sequence, but terminator is the rev_compl "))
terminator_end=input_term_end-1
terminator=Seq(str(input("Input terminator sequence (beware of lower-case carachters): ")))#Seq("AAAAGCTCTGACGGCTCCTTTTTGGAGCCGTTATCTTTTTCT".upper())
is_revcomp=str(input("is this + or - strand? [input: + or -] : "))
if is_revcomp=="+":
                    terminator=terminator.reverse_complement()
print("terminator length: ", len(terminator))
#TODO sono arrivato fino a qui ad aggiungere cose per renderlo lanciabile semi-automatico
terminator_start=terminator_end+len(terminator) #occhio che qui start e stop si riferisce al rev.compl, mentre nel antirepeat si riferisce alla direzione base


# In[32]:


# Open genome and extract the whole sequence
# find this inside the genome, and link it to the thing found by ARNold:
for record in SeqIO.parse("/shares/CIBIO-Storage/CM/scratch/tmp_projects/signorini_cas/selected_loci/"+seqid+"/"+genomename+".fa", "fasta"):
    if record.id==contigname:
        contigseq=record.seq

#check antirepeat
print("test repeat: ", antirepeat, len(antirepeat), contigseq[matched_rep_start-1:matched_rep_stop-2].reverse_complement()==antirepeat) # -1 why? Indexing diverso evidentemente col blastout, vbb
print(contigseq[matched_rep_start-1:matched_rep_stop-2].reverse_complement())
# check terminator
print("test terminator: " ,terminator, len(terminator), contigseq[terminator_end:terminator_start].reverse_complement()==terminator)


pre_tracrRNA=contigseq[terminator_end:matched_rep_stop].reverse_complement()
print("pre tracrRNA:", pre_tracrRNA, len(pre_tracrRNA))

# add mario
mario=contigseq[terminator_start:matched_rep_start-1].reverse_complement()
print("mario sequence:",  mario, len(mario))  # ho fatto i check, il -1 sopra è giusto


# ora ha una lunghezza simile, e finisce in modo simile a quello di _S. pyogenes_  TTTTTTTaltre

# In[47]:


# Now, align with the rest of the stuff to see what's up
# first up, with the repeat, lets take the rpt3, which has only got 1 mismatch

alignments = pairwise2.align.localxx(rpt3, antirepeat.reverse_complement())
print(pairwise2.format_alignment(*alignments[0]))
print("rpt length:", len(rpt3))


# In[51]:


# align antirepat + mario with repeat
alignments = pairwise2.align.globalxx( rpt3,(antirepeat+mario[:4]).reverse_complement())
print(pairwise2.format_alignment(*alignments[0]))


# In[ ]:


#no non fa nulla quel CTTT in più


# In[54]:


# I would probably say that the tracr is just, simply la parte di pre_tracrRNA la cui antirepeat si lega alla repeat, 
# che, in questo caso, è tutta.
tracrRNA = pre_tracrRNA

print("\n----------------------------------------------------------------------")
print("         tracrRNA sequence for", feature,"protein id", seqid)
print("----------------------------------------------------------------------\n")
print("- contig name:", cas_dataset[cas_dataset["Seq ID"]==seqid]["Contig"].iloc[0])
print("- Dataset: ",cas_dataset[cas_dataset["Seq ID"]==seqid]["Study"].iloc[0])
print("- Bin: ",cas_dataset[cas_dataset["Seq ID"]==seqid]["Genome Name"].iloc[0])
print("- SGB: ",cas_dataset[cas_dataset["Seq ID"]==seqid]["SGB ID"].iloc[0])
print("----------------------------------------------------------------------\n")
print("- Length: ", len(tracrRNA),"\n- position: ", terminator_end,matched_rep_stop)
print(tracrRNA)


# In[ ]:




