# Tue May 5 12:08:02 CEST 2020
# Made by L-F-S
# At the University Of Trento, Italy

import os
import subprocess
import argparse
import sys
#import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
#from Bio.SeqRecord import SeqRecord
# for clustal alignments whch i dont use anymore so TODO remove
#from Bio import AlignIO
#from Bio.Align.Applications import ClustalwCommandline
#from Bio.Align import AlignInfo
# to make logos
import logomaker as lm

# run on multiple threads
#import multiprocessing as mp
###############################################################################
#               Functions

def bulid_sequence_logos(dataset_of_flanking_sequences):
    PAMdata=pd.read_csv(outdir+dataset_of_flanking_sequences, header=None)

    # UPSTREAM Sequence mpileup and logo
    ####################################
    # Cycle over a list fo strings to write all 4 possible logos.
    counts_mat = lm.alignment_to_matrix(PAMdata[1])
    counts_mat.to_csv(outdir+"upstream.pileup")
    fig=plt.figure()
    logo=lm.Logo(counts_mat,shade_below=.5,fade_below=.5)
    logo.style_spines(visible=False)
    logo.style_spines(spines=['left', 'bottom'], visible=True)
    plt.savefig(outdir+seqid+".upstream.PAM_logo.pdf")
    plt.close()

    # Downstream
    #############################
    downstream=PAMdata[2]
    downstream_only_length=downstream[downstream.str.count("")==51] #51 because .str.count("") returns len(str)+1
    counts_mat = lm.alignment_to_matrix(downstream_only_length)
    counts_mat.to_csv(outdir+".downstream.PAM_logo.pileup")
    fig=plt.figure()
    logo=lm.Logo(counts_mat,shade_below=.5,fade_below=.5)
    logo.style_spines(visible=False)
    logo.style_spines(spines=['left', 'bottom'], visible=True)
    plt.savefig(outdir+seqid+".downstream.PAM_logo.pdf")
    plt.close()

def wrapper(seqid,feature,outdir, dataset_of_flanking_sequences):
    #subprocess.Popen("cat flanking_sequences* > dataset_flanking_sequences_of_putative_protospacers && rm flanking_sequences*", shell=True)
   # vprint("Time passed: "+str(time.time()-start_time),True,l)

    # Build sequence logos
    ##################################################################
    print("Building upstream and downstream sequence logos for cas locus",seqid)
    #TODO add input to modify length of seq logos
    bulid_sequence_logos(dataset_of_flanking_sequences)
    return







##############################################################################
if __name__=="__main__":
    main_descr="Wellcome to PAM_finder! The in-silico PAM sequence discovery essay for metagenomic data! Insert a locus id (effector id). A log file at the end of the process will be in the output folder"
    parser=argparse.ArgumentParser(description="+"*5+"\t\t"+main_descr+"\t\t"+"+"*5 )
    parser.add_argument("ID", type=str, help="sequence ID")
    parser.add_argument("-f", type=str, help="effector Cas name (default= Cas9)"\
                        , default="Cas9")
    parser.add_argument("-c", type=str, help="cas_dataset position, default=\
                        /shares/CIBIO-Storage/CM/news/users/lorenzo.signorini/5caslocitable"\
                        ,default="/shares/CIBIO-Storage/CM/news/users/lorenzo.signorini/5caslocitable/")
    parser.add_argument("-o", type=str, help="output directory, default =\
                        /shares/CIBIO-Storage/CM/news/users/lorenzo.signorini\
                        /8pamsearch/out/", default="/shares/CIBIO-Storage/CM/news/users/lorenzo.signorini/8pamsearch/out/")
    args=parser.parse_args()
    seqid =args.ID
    feature=args.f
    outdir=args.o+seqid+"/"

    if not os.path.exists(outdir):
            os.makedirs(outdir)
    dataset="dataset_flanking_sequences_of_putative_protospacers"
    wrapper(seqid,feature,outdir,dataset)
