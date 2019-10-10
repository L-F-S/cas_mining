# 23-09-2019
# Made by L-F-S
# At the University Of Trento, Italy
# see 23-09-2019_search_for_PAMS.ipynb
# USAGE:
    
    # blastsearch_for_PAM_agains_blabla <featre> <featureid>"
# Returns a .blastnlauncher file, containing several thousand blast commands to run blastn against all databases

import os
import sys 
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt 
from Bio import SeqIO
sys.path.insert(0, '/home/lorenzo.signorini/cas_mining/utils/')
import filename_discrepancies

feature=sys.argv[1]
seq_id=sys.argv[2]
datadir="/shares/CIBIO-Storage/CM/scratch/tmp_projects/signorini_cas/"

# take into account filename discrepancies
def fname_discr(dataset, cose, altrecose, probabilmente messe in dizionario)
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
    return genomename, datasetname probabilmente messe in dizionario


# Import data table
cas_dataset=pd.read_csv("/shares/CIBIO-Storage/CM/scratch/tmp_projects/signorini_cas/5caslocitable/known_"+feature+"_variants_table.csv", index_col=0)

# Extract your feature info
dataset=cas_dataset[cas_dataset["Seq ID"]==seq_id]["Study"].iloc[0]
genomename=cas_dataset[cas_dataset["Seq ID"]==seq_id]["Genome Name"].iloc[0]
SGB=cas_dataset[cas_dataset["Seq ID"]==seq_id]["SGB ID"].iloc[0]
print("\n-----------------------------------------------------\n"+feature+" id:\t"+seq_id+"\nGenome Name:\t",genomename,"\n-----------------------------------------------------\nCRISPR spacer sequences:")

######################################################################
# take into account filename discrepancies #TODO messo nella funzione sopra, che deve funzionare anche x dataset, ssample, e no bin (per il passaggio sotto
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
########################################################################

# 1 get CRISPR spacer sequence
#spacers_list_of_contig={} initial dictionary verision with contig, but I don't need the contigname

spacers_list=[]
mincedCRISPRfilename=datadir+"1crisprsearch/out/"+s3_dataset+"/"+s3_genomename+".fa.minced.out"
#print(mincedCRISPRfilename)
f=open(mincedCRISPRfilename, "r")
for line in f.readlines():
    if line.startswith("Sequence"):
        contig=line.split(" ")[1].strip("\'")
       
    if line[0]== "1" or line[0]== "2" or line[0]== "3" or line[0]== "4" or line[0]== "5" or \
    line[0]== "6" or line[0]== "7" or line[0]== "8" or line[0]== "9" or line[0]== "0":
      #  print(line)
        if not line.split("\t")[3]=="\n":
            spacers_list.append(line.split("\t")[3])
f.close()
print(spacers_list,"\n")

#2. build temporary query file for blastn search

blast_folder="/shares/CIBIO-Storage/CM/scratch/tmp_projects/signorini_cas/4exploratoryanalyses/pamsearch/out"
os.chdir(blast_folder)
tempfile=open("temp_spacer_seq", "w")
tempfile.close()
tempfile=open("temp_spacer_seq", "a")
for n, spacer in enumerate(spacers_list):

    tempfile.write(">spacer"+str(n+1)+"|"+contig+"|"+genomename+"|"+seq_id+"\n"+spacer+"\n") 
tempfile.close()
    
#3. Create a file containing the script to run a blastn search for the selected seqid against all the blastdb (one for every sample in the datasetss)
#TODO megahit mannaggia al cuore immacolato di cristo
#TODO create the dblist
def get_dbs():
    f=open("/shares/CIBIO-Storage/CM/scratch/tmp_projects/signorini_cas/list_of_samples3_name","r")
    for line in f.readlines():
        completes3name= line.strip("\n")
        s3_dataset = completes3name.split("__")[0]
        s3_sample = completes3name.split("__")[0]

        wdataset = qlcs(s3_dataset) #TODO
        megahitsample =qlcsaltro(s3_sample)

    f.close()
    return db_dictionary


db_list= {"SmitsSA_2017":["TZ_10742","TZ_99300"],"AsnicarF_2017":["MV_FEM3_t1Q14","MV_FEI5_t3Q15"]}#get_dbs()
blast_commands_for_this_sample= "# blastn commands to run Cas9 "+ seq_id + " from s3genome " + genomename+"\n"
for wdataset in db_list.keys():
    for wsample in db_list[wdataset]:
        db = "/shares/CIBIO-Storage/CM/scratch/tmp_projects/signorini_cas/4exploratoryanalyses/pamsearch/"+wdataset+"/"+wdataset+"__"+wsample+".contigs_filtered.fasta"
        blastoutfile="/shares/CIBIO-Storage/CM/scratch/tmp_projects/signorini_cas/4exploratoryanalyses/pamsearch/out/"+str(SGB)+"__"+seq_id+"__"+genomename+"__"+feature+".blastout."+wdataset+"__"+wsample
        queryfile="/shares/CIBIO-Storage/CM/scratch/tmp_projects/signorini_cas/4exploratoryanalyses/pamsearch/out/temp_spacer_seq"
        blastn_command = "blastn -out "+blastoutfile+" -outfmt \"6  qseqid sseqid pident qlen length mismatch gapopen qseq sseq sstart send evalue sstrand\" -query temp_spacer_seq -db "+db+" -evalue 0.001 -word_size 12\n"
        blast_commands_for_this_sample+=blastn_command
f=open("/shares/CIBIO-Storage/CM/scratch/tmp_projects/signorini_cas/4exploratoryanalyses/pamsearch/blastncommands/"+seq_id+"__"+genomename+"."+feature+".blastnlauncher", "w")
f.write(blast_commands_for_this_sample)
f.close()

