# Fri Apr 17 11:04:30 CEST 2020
# Made by L-F-S
# At the University Of Trento, Italy
"""
USAGE:
    originalpaths BINNAME
"""
import sys
import pandas as pd
sys.path.insert(0, '/home/lorenzo.signorini/cas_mining/utils/')
import filename_discrepancies


def print_path(bin_name, metaref=False):
    """bin (genome) names should be S3 genome names
    add optiona larg metaref=True, to print metaref path"""
    if metaref==True:
        metarefpath="/shares/CIBIO-Storage/CM/scratch/databases/MetaRefSGB/"
        datamap=pd.read_csv(metarefpath+"releases/Jan19/datamap.txt", sep="\t")
        rel_bin_path=datamap[datamap["# file_name"]==bin_name]["file_path"].iloc[0].lstrip("/")
        random_folder=rel_bin_path.split("/")[3]
        rel_anno_path=rel_bin_path.split("/")[0]+"/"+rel_bin_path.split("/")[1]+"/annotations/prokka-1.12/"+random_folder+"/"+bin_name+"/"
        return metarefpath+rel_bin_path, metarefpath+rel_anno_path
    ####
    dataset=bin_name.split("__")[0]
    bin_path="/shares/CIBIO-Storage/CM/scratch/users/e.pasolli/projects/binning/genomes_comp50_cont05/"+dataset+\
        "/"+bin_name+".fa"

    anno_path="/shares/CIBIO-Storage/CM/scratch/users/e.pasolli/projects/binning/genomes_comp50_cont05/"+dataset+\
        "/prokka/"+bin_name+"/"
    return bin_path, anno_path


if __name__=="__main__":
   print( print_path(sys.argv[1], True))
   print( print_path(sys.argv[1]))
