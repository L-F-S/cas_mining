#!/bin/bash
# Tue Jul 9 09:50:24 CEST 2019
# Made by L-F-S
# At the University Of Trento, Italy
# merges the datastes' tables back into one big table.

import os
import sys
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

sys.path.insert(0, '/home/lorenzo.signorini/utils/')
import filename_discrepancies

finaloutdir=sys.argv[1]  #5caslocitable in qst caso
feature=sys.argv[2]
outdir="/shares/CIBIO-Storage/CM/scratch/tmp_projects/signorini_cas/"+finaloutdir #TODO da cambire spesso?
CRISPRdir="/shares/CIBIO-Storage/CM/scratch/tmp_projects/signorini_cas/1crisprsearch/out"
annodir="/shares/CIBIO-Storage/CM/scratch/tmp_projects/signorini_cas/2casanno/crisprcasanno"


first=True
for dataset in ["AsnicarF_2017","BackhedF_2015","Bengtsson-PalmeJ_2015","BritoIL_2016","CM_cf","CM_madagascar","CM_periimplantitis","Castro-NallarE_2015","ChengpingW_2017","ChngKR_2016","CosteaPI_2017","DavidLA_2015","FengQ_2015","FerrettiP_2018","GeversD_2014","HMP_2012","HanniganGD_2017","HeQ_2017","IjazUZ_2017","KarlssonFH_2013","KosticAD_2015","LeChatelierE_2013","LiJ_2014","LiJ_2017","LiSS_2016","LiuW_2016","LomanNJ_2013","LoombaR_2017","LouisS_2016","NielsenHB_2014","Obregon-TitoAJ_2015","OhJ_2014","OlmMR_2017","QinJ_2012","QinN_2014","RampelliS_2015","RaymondF_2016","SchirmerM_2016","SmitsSA_2017","VatanenT_2016","VincentC_2016","VogtmannE_2016","WenC_2017","XieH_2016","YuJ_2015","ZeeviD_2015","ZellerG_2014"]:
    datasetdataframe=pd.read_csv(outdir+"/"+dataset+"/"+dataset+"_known_"+feature+"_variants_table.csv", index_col=0) #TODO da cambiare spesso?
    print(dataset+"++++++++++++++++++++++++++++++++++++++++++++")
    print(datasetdataframe.columns)
    print(datasetdataframe.index)
    if first==True:
        first=False
        bigdf=datasetdataframe
    else:
        bigdf=pd.concat([bigdf,datasetdataframe], ignore_index=True)  #TODO ignore index?
print("writing to file")
bigdf.to_csv(outdir+"/known_"+feature+"_variants_table.csv")         #TODO da cambiare spesso?
print(bigdf.shape)
