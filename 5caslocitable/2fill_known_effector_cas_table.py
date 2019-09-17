#!/bin/bash
# Fri Aug 2 09:59:19 CEST 2019
# Made by L-F-S
# At the University Of Trento, Italy

import os
import sys
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from Bio import SeqIO
from Bio import SeqIO

sys.path.insert(0, '/home/lorenzo.signorini/cas_mining/utils/')
import filename_discrepancies

####################################################################################################
#                                                 functions

def get_contig_from_tabellazza(s3genome, workingdataset,cas9ID):
    tabellapathdata=workingdataset if not workingdataset.startswith("ZeeviD") else "ZeeviD_2015"
    f=open(tabledir+"/"+tabellapathdata+"/crisprcas_hits_table_"+tabellapathdata+".csv")
#    print(s3genome)
    for line in f.readlines():
         if ","+s3genome+"," in line:#this should only happen once in the whole file. we shall belive it is so.
             line=line.strip().split(",")
             prokka=line[16]
             uniref=line[17]
            # print(uniref)
             prokkaslist=prokka.split(">")
             prokkaslist.pop(0)

             # get contigname of current effectorcas   
             for cas in prokkaslist:
                 if cas9ID in cas:
                   contig=cas.split("\t")[0]
                   f.close()
                   return contig


####################################################################################################
#                                            MAIN
####################################################################################################
wdataset=sys.argv[1] #TODO switch test dataset
s3dataset=filename_discrepancies.s3(wdataset)
feature=sys.argv[2]
tabledir="/shares/CIBIO-Storage/CM/scratch/tmp_projects/signorini_cas/3tabellazza"
outdir="/shares/CIBIO-Storage/CM/scratch/tmp_projects/signorini_cas/5caslocitable/"+s3dataset
datadir="/shares/CIBIO-Storage/CM/scratch/tmp_projects/signorini_cas/" 
annodir="/shares/CIBIO-Storage/CM/scratch/tmp_projects/signorini_cas/2casanno/crisprcasanno"

DF=pd.read_csv(outdir +"/"+s3dataset+"_known_"+feature+"_variants_table.csv", index_col=0) #input e' come 3tabellazza, ma solo il subset con Cas9. piu piccolo e agile.

outdf=pd.DataFrame(columns=["Seq ID","Seq Description","Seq","Contig","Genome Name","Study","Sample Name","SGB ID"])
newindexmitocca=0

for index, genome in DF.iterrows(): 
    genomename=genome["Genome Name"] #TODO Be careful for filename discrepancies, especially with ZeeviD files and with _megahit_ underscores!
    SGB=genome["SGB ID"]
    print("\n-----------------------------------------------------\nGenome Name:",genome["Genome Name"],"\n-----------------------------------------------------\n")
    working_genomename, working_dataset=filename_discrepancies.get_originalsamplename_froms3name_of_genome(genomename, s3dataset) 

    path="/shares/CIBIO-Storage/CM/scratch/tmp_projects/epasolli_darkmatter/allcontigs/"+working_dataset+"/metabat/genomes_comp50_cont05/prokka/"+working_genomename    
    filename=path+"/"+working_genomename+".faa"
    for record in SeqIO.parse(filename, "fasta"):
        if feature in record.description:
           # print("UN CAS9!")
            # devo comunque chiaramente mettere tutto qui
            outdf.at[newindexmitocca,"Seq ID"]=str(record.id)
            outdf.at[newindexmitocca,"Seq Description"]=str(record.description)
            outdf.at[newindexmitocca,"Seq"]=str(record.seq)
            outdf.at[newindexmitocca,"Genome Name"]=genomename
            outdf.at[newindexmitocca,"Study"]=s3dataset
            outdf.at[newindexmitocca,"Sample Name"]=genome["Sample Name"]
            outdf.at[newindexmitocca,"SGB ID"]=SGB
            outdf.at[newindexmitocca,"Contig"]=get_contig_from_tabellazza(genomename,working_dataset,record.id)
            newindexmitocca+=1     
#
outdf.to_csv(outdir+"/"+s3dataset+"_known_"+feature+"_variants_table.csv")  #this is overwriting. occhio
