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

def get_SGB_info(SGB):
    f=open("/shares/CIBIO-Storage/CM/scratch/tmp_projects/signorini_cas/S4Segata.csv")
#    print(s3genome)
    for line in f.readlines():
        if ","+str(SGB)+"," in line:#this should only happen once in the whole file. we shall belive it is so.
            line=line.strip().split(",")
            uSGB=line[4]
#            print(uSGB)
            level_estim=line[5]
            estim=line[6]
            f.close()
            return uSGB, estim, level_estim


###############################################################################################################
#                           MAIN                                                                   
feature=sys.argv[1]
outdir="/shares/CIBIO-Storage/CM/scratch/tmp_projects/signorini_cas/5caslocitable/"

DF=pd.read_csv(outdir +"known_"+feature+"_variants_table.csv", index_col=0) #input e' come 3tabellazza, ma solo il subset con Cas9. piu piccolo e agile.

#outdf=pd.DataFrame(columns=["Seq ID","Seq Description","Seq","Contig","Genome Name","Study","Sample Name","SGB ID"])


for index, genome in DF.iterrows():
        genomename=genome["Genome Name"] #TODO Be careful for filename discrepancies, especially with ZeeviD files and with _megahit_ underscores!
        SGB=genome["SGB ID"]
#        print("\n-----------------------------------------------------\nGenome Name:",genome["Genome Name"],"\n-----------------------------------------------------\n")
        uSGB, estim, level_estim =get_SGB_info(SGB)
        DF.at[index,"uSGB"]=uSGB
        DF.at[index,"Level of estimated taxonomy"]=level_estim
        DF.at[index,"Estimated taxonomy"]=estim
#        print("trovati: usgb, estim, level_estim:", uSGB, estim, level_estim)
DF.to_csv(outdir+"known_"+feature+"_variants_table.csv")  #this is NOT  overwriting. occhio

