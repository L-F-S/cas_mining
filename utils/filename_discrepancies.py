#!/bin/bash
# 06/06/2019 non e' vero data sbagliata lol
# Made by L-F-S
# At the University Of Trento, Italy
# Import this into your other scripts to help you managing tricky filenames

import os
import sys
import numpy as np
import pandas as pd

def Schirmer_2016(genomename,from_anno=False):
    names={"SchirmerM_2016__G88886__bin.2":"SchirmerM_2016__G88886__bin","SchirmerM_2016__G88887__bin.10":"SchirmerM_2016__G88887__bin","SchirmerM_2016__G88911__bin.27":"SchirmerM_2016__G88911__bin","SchirmerM_2016__G89138__bin.29":"SchirmerM_2016__G89138__bin"}
    if not from_anno:
        if genomename in names.keys():
#            print("Maledetto schirmer che hai i nomi diversi!")
            return names[genomename]
        else:
            return genomename 

def analysis_has_megahit(analysis):
    """takes as input: "prokka", "pilercr", "minced", "old_anno" returns a BOOL indicating whether this has megahit in the name"""
    return True if analysis=="prokka" or analysis == "old_anno" else False

def dataset_has_megahit(dataset,sample):# in realta' e' solo zeevid (senz aa e b)
#funzionicchia, ma se sample  e' in names3 e dataset no , scazza
 # attenzione!!! OhJ_2014 ce l'ha per alcuni cosi si e per altri no!!
  # attenzione!!! OhJ_2014 ce l'ha per alcuni cosi si e per altri no!!
    """ returns True if the given sample has the megahit interfix """
    if dataset == "BengtssonPalmeJ_2015" or dataset ==  "OhJ_2014" or dataset == "IjazUZ_2017" or  dataset == "LiJ_2014" or  dataset == "Castro-Nallar" or dataset == "Castro-NallarE_2015" or  dataset == "RaymondF_2016" or  dataset == "SmitsSA_2017" or  dataset == "VatanenT_2016" or dataset == "VincentC_2016" or dataset == "WenC_2017" or dataset == "ZeeviD_2015_A" or dataset == "ZeeviD_2015_B":
        
  #      if dataset in ["IjazUZ_2017",  "Castro-NallarE_2015",  "SmitsSA_2017", "VincentC_2016", "ZeeviD_2015_A", "ZeeviD_20015_B","ZeeviD"]:
   #         return True
        megahit_name=change_to_megahit(sample)
        try:
            f=open("/shares/CIBIO-Storage/CM/scratch/tmp_projects/signorini_cas/2casanno/justminced/"+dataset+"/"+megahit_name+".crisprcas.gff.minced")
            f.close()
            return True
        except:
            return False
    return False
            
  #  else:
       # return true if "MET0008" in sample else False
def s3(dataset, r=False): 
    """returns the S3table name of the given <dataset>, if different. Otherwise, returns <dataset>.
    if r=True, does the reverse: takes as input the S3name and returns the working name"""
    different_dataset_names = {"BengtssonPalmeJ_2015" : "Bengtsson-PalmeJ_2015", "LawrenceA_2015": "DavidLA_2015","CM_caritro":"FerrettiP_2018", "ZeeviD_2015_A":"ZeeviD_2015", "ZeeviD_2015_B":"ZeeviD_2015"}
    s3_names=["Bengtsson-PalmeJ_2015",  "DavidLA_2015","FerrettiP_2018", "ZeeviD_2015", "ZeeviD_2015"]
    working_names = ["BengtssonPalmeJ_2015" , "LawrenceA_2015","CM_caritro", "ZeeviD_2015_A", "ZeeviD_2015_B"] 
    if not r:
        if not dataset in different_dataset_names.keys():
            return dataset
        else:
            return different_dataset_names[dataset]
    else:
        if not dataset in different_dataset_names.values():  #TODO does not work for zeevid, se ce zeevid te ne se i occupato nello script
            return dataset
        else:
            index=s3_names.index(dataset)
            return working_names[index]


def change_to_megahit(filename):
    """takes as input the filename in the dataset__sample__bin.n  format, and returns a string with the addition of the megahit interfix:
        dataset__sample_megahit__bin.n"""
    l= filename.split("__")
    sample_megahit=l[0]+"__"+l[1]+"_megahit__"+l[2]
    return sample_megahit 

def is_this_sample_in_ZeeviD(dataset,filename):
    """takes as input the working dataset name and the filename in the ZeeviD_2015__sample__bin.n  format, and returns it in the 
    ZeeviD__2015_LETTER__sample_megahit__bin.n format"""
    # os.chdir("/shares/CIBIO-Storage/CM/scratch/tmp_projects/signorini_cas/2casanno/justminced/ZeeviD_2015_A")
    letter=dataset[-1]
    f=open("/shares/CIBIO-Storage/CM/scratch/tmp_projects/signorini_cas/all_ZeeviD_2015_B")
    lines=f.readlines()
    f.close()
    strippedfilename=filename.split("__")[1]+"_megahit"
    strippedlines=[line.split("__")[1] for line in lines]
    if strippedfilename in strippedlines:
        return "ZeeviD_2015_B__"+strippedfilename+"__"+filename.split("__")[2], "B" 
    return "ZeeviD_2015_A__"+strippedfilename+"__"+filename.split("__")[2], "A"



def get_originalsamplename_froms3name_of_genome(genomename, dataset):
    """input: genome name and dataset name as in table S3, output: genomename and dataset name as in epasolli/darkmatter/allcontigs. """
    if dataset.startswith("ZeeviD"):
        working_genomename=change_to_megahit(genomename)
        try:
            working_dataset="ZeeviD_2015_B"
            working_genomename=working_genomename.replace("ZeeviD_2015", "ZeeviD_2015_B")
            prokkafile_ofsample=open("/shares/CIBIO-Storage/CM/scratch/tmp_projects/signorini_cas/2casanno/justminced/ZeeviD_2015_B/"+working_genomename+".crisprcas.gff.minced")
            prokkafile_ofsample.close()
        except:
            working_genomename=working_genomename.replace("ZeeviD_2015_B", "ZeeviD_2015_A")
            working_dataset="ZeeviD_2015_A"
    else:
        working_dataset=s3(dataset,r=True)
        working_genomename=genomename.replace(dataset,working_dataset)
    if dataset_has_megahit(working_dataset,working_genomename):
        working_genomename=change_to_megahit(working_genomename)
    return working_genomename, working_dataset

def master_filename_discrepancies(dataset,genomename):
    # dataset and genomename are working dataset and genomename
    if dataset.startswith("ZeeviD"):  #TODO guarda se funziona quando finisci di girare minced_runner
        s3_dataset="ZeeviD_2015"
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
        dataset=s3(s3_dataset,r=True) # get working dataset name from s3 dataset name
        genomename=genomename.replace(s3_dataset, working_dataset)

    if dataset_has_megahit(dataset,genomename):
        genomename=change_to_megahit(genomename)

    samplename = genomename.split("__")[1]
                                      
    return dataset, genomename, samplename, s3_dataset, s3_genome
#if __name__=="__main__":
#

