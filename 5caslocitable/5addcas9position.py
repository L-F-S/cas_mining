# Wed Aug 7 08:41:17 CEST 2019
# Made by L-F-S
# At the University Of Trento, Italy

import os
import sys
import numpy as np
import pandas as pd

sys.path.insert(0, '/home/lorenzo.signorini/cas_mining/utils/')
import filename_discrepancies
print("MERGE ALL DATASETS TOGETHER BEFORE RUNNING THIS STEP!!!!!!!!!!!!!!!!!!!!")
#megia iniseme all the dataset before running this step!!

#TODO : to be done after having merged all the datasets because it works on the big one, in teoria, ma in pratica fose no shalla si puo fare anche per meno

##############################################################################################################

def get_feature_pos(seqid, genomename):
    features_of_genome=hits_table.loc[hits_table["Genome Name"]==genomename, "prokka_cas"].iloc[0]
    features_of_genome=features_of_genome.strip().lstrip(">").split(">")
    for single in features_of_genome:
        if seqid in single:
            split_single_side_feature=single.split("\t")
            print(split_single_side_feature)
            idd, start, end=split_single_side_feature[8].split(";")[0], split_single_side_feature[3],split_single_side_feature[4]
    tupla=(idd,start,end)
    return tupla 

###############################################################################################################
#                           MAIN                                                                   
feature=sys.argv[1]
outdir="/shares/CIBIO-Storage/CM/scratch/tmp_projects/signorini_cas/5caslocitable/"
tabledir="/shares/CIBIO-Storage/CM/scratch/tmp_projects/signorini_cas/3tabellazza"
DF=pd.read_csv(outdir +"known_"+feature+"_variants_table.csv", index_col=0) #input e' come 3tabellazza, ma solo il subset con Cas9. piu piccolo e agile.
hits_table=pd.read_csv(tabledir+"/crisprcas_hits_table.csv", index_col=0)


previous_dataset="pippo"
prokka_cas9_list=[]
for index, genome in DF.iterrows():
        genomename=genome["Genome Name"] #TODO Be careful for filename discrepancies, especially with ZeeviD files and with _megahit_ underscores!
        if previous_dataset!=genome["Study"]:
            print(genome["Study"])
            previous_dataset=genome["Study"]
        prokka_cas9 = get_feature_pos(genome["Seq ID"], genomename)
        if np.shape(prokka_cas9)!=(3,):
            print(prokka_cas9,"\n\n------------------")
        prokka_cas9_list.append( prokka_cas9)
#        print("trovati: usgb, estim, level_estim:", uSGB, estim, level_estim)
DF["prokka_cas9"]=pd.Series(prokka_cas9_list, index=DF.index)
DF.to_csv(outdir+"known_"+feature+"prokkacas9_variants_table.csv")  #this is NOT  overwriting. occhio

