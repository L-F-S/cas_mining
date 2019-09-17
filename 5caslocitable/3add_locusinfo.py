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

#TODO big todo: pilercr has missing 'stop' position in handmade gff annotation :,(.
# USAGE: python 3add_blbabla.py dataset feature [Cas9, Cpf1..]
####################################################################################################
#                                                 functions

def get_feature_pos(side_feature, contig, s3genome):
    #HERE WE ARE INSIDE ONE SINGLE ROW (ONE GENOME)
    # and we are considering one side  feature( es: CRISPR) for one main feature (es: Cas9) on one contig

    #returns a series of tuples. Every tuple is in the shape (side featureID, side feature start, side feature end)
    # side features are CRISPR arrays or non-effector Cas proteins
    # check name nel caso particoalre di estrarre piu cas dalla riga dei cas:
    column_name=side_feature
    if side_feature.startswith("prokka_cas"):
        column_name="prokka_cas"

    features_of_genome=hits_table.loc[hits_table["Genome Name"]==s3genome, column_name].iloc[0]#.iloc[0] to extract the first (only!) element of that pd.Series!

    series_of_tuples=[]

    if type(features_of_genome)== str:
        features_of_genome=features_of_genome.strip().lstrip(">").split(">")  


        for single_sidefeature in features_of_genome:
            if contig in single_sidefeature:
                if not "Cas9" in single_sidefeature:
                    split_single_side_feature=single_sidefeature.split("\t")  
#                    print(split_single_side_feature, len(split_single_side_feature))
                    if side_feature.startswith("prokka"):
                        if side_feature.startswith("prokka_cas1"):   #handle cas1 and cas2
                            if "Cas1" in single_sidefeature:
                                idd, start, end=split_single_side_feature[8].split(";")[0], split_single_side_feature[3],split_single_side_feature[4]
                                series_of_tuples.append((idd, start, end))
                        if side_feature.startswith("prokka_cas2"):
                            if "Cas2" in single_sidefeature:
                                idd, start, end=split_single_side_feature[8].split(";")[0], split_single_side_feature[3],split_single_side_feature[4]
                                series_of_tuples.append((idd, start, end))
    
                    else:                        
                        idd, start, end=split_single_side_feature[7].split(";")[0], split_single_side_feature[3], split_single_side_feature[4]
                        series_of_tuples.append((idd, start, end))
                
                    

        return series_of_tuples if len(series_of_tuples)>0 else None
    else: #simply no side_feature for that genome.
        return None 


def dataframe_row_iterator():
    # not passing any input because the main is defined globally.
    # iterate over lines of cas effector
    # return several Series of tuples of ID, start and end position of given feature. if feature is not present, returns nan
    list_of_series_of_features={}
    for index, cas9series in DF.iterrows(): 
        genomename=cas9series["Genome Name"] #TODO Be careful for filename discrepancies, especially with ZeeviD files and with _megahit_ underscores!
        SGB=cas9series["SGB ID"]
        contig=cas9series["Contig"]
        print("\n-----------------------------------------------------\nGenome Name:",cas9series["Genome Name"],"\n-----------------------------------------------------\n")
        working_genomename, working_dataset=filename_discrepancies.get_originalsamplename_froms3name_of_genome(genomename, s3dataset) 

        for side_feature in ["pilercr_CRISPR","minced_CRISPR","prokka_cas1", "prokka_cas2"]: #TODO cambia qui la feature che vuoi trovare!
            print(side_feature)
            if not side_feature in list_of_series_of_features.keys():
                list_of_series_of_features[side_feature]=[]
            
            # get[ (ID, start, end), (ID, start, end)..] for the side feature(s) of that Cas9 row!
            list_of_series_of_features[side_feature].append(get_feature_pos(side_feature, contig, genomename))
 #       print( list_of_series_of_features)
        
    # TRANSFORM INTO SERIES:
    for side_feature in list_of_series_of_features.keys():
        list_of_series_of_features[side_feature]=pd.Series(list_of_series_of_features[side_feature], index=DF.index)
    return list_of_series_of_features 


####################################################################################################
#                                            MAIN
####################################################################################################
# stavolta il wrapper che itera le righe non e' nella funzione ma e' definito globale nel main. 
#così vedi bene che stai facendo lo stesso ciclo più volte

wdataset=sys.argv[1] 
s3dataset=filename_discrepancies.s3(wdataset)
feature=sys.argv[2]  # es: Cas9, Cpf1..
tabledir="/shares/CIBIO-Storage/CM/scratch/tmp_projects/signorini_cas/3tabellazza"
outdir="/shares/CIBIO-Storage/CM/scratch/tmp_projects/signorini_cas/5caslocitable/"+s3dataset
datadir="/shares/CIBIO-Storage/CM/scratch/tmp_projects/signorini_cas/" 
annodir="/shares/CIBIO-Storage/CM/scratch/tmp_projects/signorini_cas/2casanno/crisprcasanno"

DF=pd.read_csv(outdir +"/"+s3dataset+"_known_"+feature+"_variants_table.csv", index_col=0) 
hits_table=pd.read_csv(tabledir+"/crisprcas_hits_table.csv", index_col=0)

print(DF.columns, DF.shape)

#. add anything you like:
dictionary_of_series_of_sidefeatures = dataframe_row_iterator()  #TODO remember to change passed output according to what you want.
for side_feature in ["pilercr_CRISPR","minced_CRISPR","prokka_cas1","prokka_cas2"]:
    DF[side_feature]=dictionary_of_series_of_sidefeatures[side_feature]
#TODO devi discernerle ancora le prokkacas




# print stuff and write to file
print(DF.columns, DF.shape)


print(DF.head())

DF.to_csv(outdir+"/"+s3dataset+"_known_"+feature+"_variants_table.csv")  #this is overwriting. occhio
