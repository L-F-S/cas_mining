# Fri Dec 6 10:18:34 CET 2019
# Made by L-F-S
# At the University Of Trento, Italy


import os
import sys
import pandas as pd
import argparse

import locus

def get_ID_info(seqid, feature,v):
    outdir="/shares/CIBIO-Storage/CM/scratch/tmp_projects/signorini_cas/Nstep" #TODO add nstepdir
    datadir="/shares/CIBIO-Storage/CM/scratch/tmp_projects/signorini_cas/" #TODO add nstepdir
    annodir="/shares/CIBIO-Storage/CM/scratch/tmp_projects/signorini_cas/2casanno/crisprcasanno"
    caslocus=locus.locus(seqid,feature)
    cas_dataset=pd.read_csv("/shares/CIBIO-Storage/CM/scratch/tmp_projects/signorini_cas/5caslocitable/known_"+feature+"prokkacas9_variants_table.csv", index_col=0)
    caslocus.set_SGB(cas_dataset)
    caslocus.set_contigname(cas_dataset)
    caslocus.set_genomename(cas_dataset)
    caslocus.set_datasetname(cas_dataset)
    caslocus.set_sequence(cas_dataset)
    print("-"*80+"\n"," "*int((80- len("Sequence ID"))/2)," Sequence ID",\
          " "*int((80- len("Sequence ID"))/2)+"\n",\
          " "*int((80- len("Sequence ID"))/2)+ seqid+\
          " "*int((80- len("Sequence ID"))/2)+"\n"+"-"*80+"\n")
    sequence_length=len(caslocus.seq)
    print("-> Contig Name:\t",caslocus.contigname)
    print("-> Genome Name:\t",caslocus.genomename)
    print("-> Dataset Name:\t",caslocus.datasetname)
    print("-> SGB ID:\t",caslocus.SGB)
    print("-> Sequence length:\t",sequence_length)
    tmp=locus.CRISPRarray(feature=feature, contigname=caslocus.contigname,genomename=caslocus.genomename,datasetname=caslocus.datasetname)
    tmp.get_CRISPR_array(v=v)
    caslocus.CRISPRarray=tmp #Magari slightly redundant to have that locus
    caslocus.fetch_positions(cas_dataset)
    print("-> Positions:\t",caslocus.positions)
    print("TODO: Finding tracr_RNA")
    caslocus.tracrRNA="AAATGTGCATCGTACGTCAGCTGATCGTGCTACGTACGATCGATCG"
    print("\n----------------------------------------------------------------------")
    print("         tracrRNA sequence for", feature,"protein id", seqid)
    print("----------------------------------------------------------------------\n")
    print("- Length: ", len(caslocus.tracrRNA),"\n- position: ", 1,2)
    print(caslocus.tracrRNA)
    #TODO add the written output


if __name__=="__main__":
    parser=argparse.ArgumentParser(description="Print information about a specific effector Cas Sequence ID. (Input: <Seq ID>)")
    parser.add_argument("-v", action="store_true", help="verbose output")
    parser.add_argument("ID", type=str, help="sequence ID")
    parser.add_argument("feature", type=str, help="effector Cas name (i.e. Cas9)")
    args=parser.parse_args()

    seqid =args.ID
    feature=args.feature

    get_ID_info(seqid, feature,args.v)
