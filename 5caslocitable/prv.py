#!/bin/bash
# 01/07/2019
# Made by L-F-S
# At the University Of Trento, Italy
# Takes as input $data/S3Segata/S3Segata_<dataset>.csv, (which are the tables from the paper, split by datasets) and returns a new table with 4 more columns: CRISPR_pilercr, CRISPR_minced, CRISPR_prokka, cas_prokka, containing gff3 formatted annotations of interesting CRISPRCas  things. Saves the new tables in $data/3tabellazza/<dataset>/cripr_hits_table_<dataset>.csv
#
#       USAGE:
#       
#       python tabellazza.py <datasetname>
#

import os
import sys
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import time
sys.path.insert(0, '/home/lorenzo.signorini/utils/')
import filename_discrepancies
overall_time = time.time()

# temp returns a dictionary which has the column names as keys and a pd.Series object of length df.shape[0] and with index=df.index:
#    return {key:value for (key,value) in [(i+1,pd.Series(np.random.randn(df.shape[0]), index=df.index)) for i in range(np.random.randint(2,6))]


def gather_cas_info_from_uniref(df,dataset,S3_alternative_dataset_name):
    start_time=time.time()
    print("adding uniref cas entries...")
    unirefannodir="/shares/CIBIO-Storage/CM/scratch/tmp_projects/signorini_cas/6mainrealtaeprimadel3unirefanno"
    cas_from_uniref_for_all_samples=[]
    for i, row in df.iterrows():
        sample=row["Genome Name"]
        #### take into account megahitmegahit  filename discrepancies####
        # nb: visto che df è il crisprcas_hits_table_<dataset> gia processato (xke qst e un aggiunta posteriore, il genome name al suo interno è gia stato cambiato in <s3_datasetname>
        # quindi devo aggiungere questo if per checkare il megahit nei dataset merdosi  tipo bengtsonpalmer che hanno il s3 name diverso ,e hanno megahit dentro alcuni sample edentro altri no
        # quindi devo aggiungere questo if per checkare il megahit nei dataset merdosi  tipo bengtsonpalmer che hanno il s3 name diverso ,e hanno megahit dentro alcuni sample edentro altri no
        samplenameformegahit=sample
        if dataset != S3_alternative_dataset_name:
            samplenameformegahit=sample.replace(S3_alternative_dataset_name,dataset)
        if filename_discrepancies.dataset_has_megahit(dataset, samplenameformegahit):
            sample=filename_discrepancies.change_to_megahit(sample)
        if sample.startswith("SchirmerM_2016"):
            sample=filename_discrepancies.Schirmer_2016(sample)
        else:
#            print("SONO QUA VA BENE!", os.getcwd(), sample)
            sample=sample
        #################################################

        unireffile_ofsample=open(unirefannodir+"/"+S3_alternative_dataset_name+"/"+sample+".annotated.unirefanno")
        line=unireffile_ofsample.readline()
        unireffile_ofsample.close()
        cas_from_uniref_for_all_samples.append(line.rstrip("\n"))
    print("Done. Elapsed time: ", str(time.time()-start_time))
    return pd.Series(cas_from_uniref_for_all_samples, index=df.index)

##################################################################################################################
#
#                                               MAIN        
#
#################################################################################################################

def main(from_minced):
    # dataset means working dataset
    dataset=sys.argv[1] 
    print("++++++++++++++++\nadding CRISPR cas info to dataset: ", dataset)
    outdir="/shares/CIBIO-Storage/CM/scratch/tmp_projects/signorini_cas/3tabellazza/"+dataset 
    CRISPRdir="/shares/CIBIO-Storage/CM/scratch/tmp_projects/signorini_cas/1crisprsearch/out"
    annodir="/shares/CIBIO-Storage/CM/scratch/tmp_projects/signorini_cas/2casanno"
    S3_alternative_dataset_name=filename_discrepancies.s3(dataset)
    input_table="/shares/CIBIO-Storage/CM/scratch/tmp_projects/signorini_cas/3tabellazza/" + dataset +"/crisprcas_hits_table_"+dataset+".csv"    # in queste tables, il datasetname del nome della tabella è il working, quello dei sample è invece il s3, tranne per Zeevid, che è . leva  tutti i punti dall uno al 4 a ggiungi un 5 punto
    print("S3 datasetname = ", S3_alternative_dataset_name)
    df=pd.read_csv(input_table, index_col=0) 
    print(df.columns)
    print("initial dataset (rows, cols): ", df.shape)
    print("overall S3 index of current dataset's bins: from ",df.index[0],"to: ", df.index[-1])
   
   # add column. 
    unirefcas=gather_cas_info_from_uniref(df,dataset,S3_alternative_dataset_name)
    print(unirefcas.shape, " entries added")
    df["uniref_cas"]=unirefcas


# Print summary info and write to file 
#    print("--------------------\nN of bins doublecheck:\npilercr:\t{} bins\nminced: \t{} bins\nprokkaCRISPR:\t{} bins\nprokkaCas:\t{} bins\n-------------------".format(pilercr_bins,minced_bins,prokka_crispr_bins,prokka_cas_bins)) 
 #   if not pilercr_bins==minced_bins==prokka_cas_bins==prokka_crispr_bins:
  #      print("ERROOOORRRR!!!!! N of bins varying across algorithms")
    print(df.columns)
    print("new dataset (rows, cols): ", df.shape)
    name="crisprcas_hits_table_"+dataset+".csv"
    print("writing to file: " + name +"\n++++++++++++++++")
    if dataset.startswith("ZeeviD"):
        df.to_csv("/shares/CIBIO-Storage/CM/scratch/tmp_projects/signorini_cas/3tabellazza/ZeeviD_2015/crisprcas_hits_table_ZeeviD_2015.csv", quoting=3, escapechar="\\")
    else:
        df.to_csv(outdir+"/"+name, quoting=3, escapechar="\\")
    print("DONE. Total elapsed time: ", time.time()-overall_time)
    return

main(from_minced=True)
