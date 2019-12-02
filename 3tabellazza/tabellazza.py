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
sys.path.insert(0, '/home/lorenzo.signorini/cas_mining/utils/')
import filename_discrepancies
overall_time = time.time()

# temp returns a dictionary which has the column names as keys and a pd.Series object of length df.shape[0] and with index=df.index:
#    return {key:value for (key,value) in [(i+1,pd.Series(np.random.randn(df.shape[0]), index=df.index)) for i in range(np.random.randint(2,6))]

def make_a_gff_like_annotation_from_pilercr_output(pilercrfile):
    # ma forse anche un po' piu accurato del minced, tipo +/- per dire l'orientazione? ma sticazzi in
    # realta' .
    # minced entry:
    # NODE_152_length_7120_cov_13.0000_ID_1969        minced:0.4.0   repeat_region   56      683     10      .     .ID=CRISPR1;rpt_type=direct;rpt_family=CRISPR;rpt_unit_seq=ATTGTAGTTCCCTAATTTTTTTTGGTATGTTATAAT
    siamo_al_summary_by_position=False
    siamo_alla_line_che_ci_interessa=False
    stiamo_guardando_le_lines_di_un_node=False
    n=0
    gff_line_totale_che_comprende_piu_crispr_se_ci_sono = ''
    # GFF3 file format has the follwoing nine entries
    gff_sequid='.'
    gff_source="pilercr1.06"#TODO pilercr:1.06 ?
    gff_type='repeat_region'
    gff_start='.' #= Position entry in .pilecr.out file
    gff_end='.'  # position+length of .pilercr.ou file
    gff_score='.' # nel nostro caso è il numero di copie # copies
    gff_strand='.'
    gff_phase='.'   # not interesting for us. It shall remain empty.
    gff_attributes = '.' #  list of feature attributes in the format tag=value;tag2=value2
    ######################################################################################
    # the following are tags inside the attributes column:
    attr_ID="ID="#TODO puo esee anche2 (0 3?)
    attr_rpt_type="rpt_type=direct"
    attr_rpt_family="rpt_family=CRISPR"
    attr_rpt_unit_seq="rpt_unit_seq="
    for line in pilercrfile.readlines():
        if line.startswith("SUMMARY BY POSITION"):
            siamo_al_summary_by_position=True
        if siamo_al_summary_by_position == True:
            if line.startswith(">"):
                n+=1 #to count the number of CRISPR inside the file
                stiamo_guardando_le_lines_di_un_node=True
                gff_sequid=line.lstrip(">").rstrip("\n")
                attr_ID="ID=CRISPR"+str(n)
            if siamo_alla_line_che_ci_interessa:
    #            print(line.split())
                gff_start=line.split()[2]
                gff_end=str(int(gff_start)+int(line.split()[3])-1)
                gff_score=line.split()[4]
               # gff_strand=line.split()[7]
                attr_rpt_unit_seq+=line.split()[7]
                gff_line_totale_che_comprende_piu_crispr_se_ci_sono +=  ">"+gff_sequid+"\t"+gff_source+"\t"+gff_start+"\t"+gff_end+"\t"+gff_score+"\t"+gff_strand+"\t"+gff_phase+"\t"+attr_ID+";"+attr_rpt_type+";"+attr_rpt_family+";"+attr_rpt_unit_seq+"\t"
   #             print(gff_line_totale_che_comprende_piu_crispr_se_ci_sono)
                siamo_alla_line_che_ci_interessa=False

            if line.startswith("="):
                siamo_alla_line_che_ci_interessa=True

    return gff_line_totale_che_comprende_piu_crispr_se_ci_sono

def gather_crispr_info_from_PilerCR(df,dataset,CRISPRdir,S3_alternative_dataset_name):
    # since tehre is no gff formatted annotation for pilercr, we must do it ourselves. Let us copy the minced format.
    start_time=time.time()
    print("adding pilercr CRISPR entries...")
    os.chdir(os.path.join(CRISPRdir,S3_alternative_dataset_name))
    CRISPR_from_pilercr_for_all_samples=[]
    count_bins=0
    for i, row in df.iterrows():
        count_bins+=1
        sample_s3name=row["Genome Name"] #rather sloppy of me to call this sample, tecnically it is a bin
   #    sample=sample_s3name
   #     if dataset != S3_alternative_dataset_name:
   #         sampletemp=sample.lstrip(S3_alternative_dataset_name)
   #         sample=dataset+"__"+sampletemp
        pilercrfile_ofsample=open(CRISPRdir+"/"+S3_alternative_dataset_name+"/"+sample_s3name+".fa.pilercr.out")
        CRISPR_from_pilercr=make_a_gff_like_annotation_from_pilercr_output(pilercrfile_ofsample)
        pilercrfile_ofsample.close()
        CRISPR_from_pilercr_for_all_samples.append(CRISPR_from_pilercr)
    print("Done. Elapsed time: ", str(time.time()-start_time))
    return pd.Series(CRISPR_from_pilercr_for_all_samples, index=df.index), count_bins

def gather_crispr_info_from_minced(df,dataset, CRISPRdir,S3_alternative_dataset_name): # WORKNIG! ma, TODO: occhio agli input di ZeeviD e agli output dei cosi sbajad
    start_time=time.time()
    print("adding minced CRISPR entries...")
    os.chdir(os.path.join(CRISPRdir,S3_alternative_dataset_name))
    CRISPR_from_minced_for_all_samples=[]
    count_bins=0
    for i, row in df.iterrows():
        count_bins+=1
        sample=row["Genome Name"]
   #     if dataset != S3_alternative_dataset_name:
   #         sampletemp=sample.lstrip(S3_alternative_dataset_name)
   #         sample=dataset+"__"+sampletemp
        mincedfile_ofsample=open(CRISPRdir+"/"+S3_alternative_dataset_name+"/"+sample+".fa.minced.out.gff")
        CRISPR_from_minced=""
        for line in mincedfile_ofsample.readlines():
            CRISPR_from_minced +=">"+line.rstrip("\n")
            if "," in line:
                print("ABOOOOOOORT!!")
        mincedfile_ofsample.close()
        CRISPR_from_minced_for_all_samples.append(CRISPR_from_minced)
    print("Done. Elapsed time: ", str(time.time()-start_time))
    return pd.Series(CRISPR_from_minced_for_all_samples, index=df.index), count_bins

def gather_crispr_info_from_prokka(df,annodir,dataset,S3_alternative_dataset_name,is_ZeeviD_A,is_ZeeviD_B): # WORKNIG! ma, TODO: da fare TODO hannno megahit nel nome
    start_time=time.time()
    print("adding prokka CRISPR entries...")
    os.chdir(annodir+"/justminced/"+dataset)
    CRISPR_from_prokka_for_all_samples=[]
    count_bins=0
    for i, row in df.iterrows():
        count_bins+=1
        sample=row["Genome Name"]
        if dataset.startswith("ZeeviD"): #nomenclature exception
            os.chdir("../ZeeviD_2015_B")
            sample=filename_discrepancies.change_to_megahit(sample)
            try:
                sample=sample.replace("ZeeviD_2015", "ZeeviD_2015_B")
                prokkafile_ofsample=open(sample+".crisprcas.gff.minced")
               # print(sample, "uweuwe")
            except:
                os.chdir("../ZeeviD_2015_A")
                sample=sample.replace("ZeeviD_2015_B", "ZeeviD_2015_A")
                prokkafile_ofsample=open(sample+".crisprcas.gff.minced")
               # print(sample, "uweuwe")

#                sample, letter=filename_discrepancies.is_this_sample_in_ZeeviD(dataset,sample)
 #               if dataset[-1] != letter:
  #                  continue
        else:
            if dataset != S3_alternative_dataset_name:
                sample=sample.replace(S3_alternative_dataset_name, dataset)
            if filename_discrepancies.analysis_has_megahit("prokka"):    # True always. redundant as fuck.
                if filename_discrepancies.dataset_has_megahit(dataset,sample):  # False per LawrenceA_2015
                    sample_megahit=filename_discrepancies.change_to_megahit(sample)
                    prokkafile_ofsample=open(sample_megahit+".crisprcas.gff.minced")
                else:                                                    # True per LawrenceA_2015
                    prokkafile_ofsample=open(sample+".crisprcas.gff.minced")
            else:
                print("This shouldn't even happen")
                break
        CRISPR_from_prokka=""
        for line in prokkafile_ofsample.readlines():
            CRISPR_from_prokka +=">"+line.rstrip("\n")
            if "," in line:
                print("ABOOOOOOORT!!")
        prokkafile_ofsample.close()
        CRISPR_from_prokka_for_all_samples.append(CRISPR_from_prokka)
    print("Done. Elapsed time: ", str(time.time()-start_time))
    return pd.Series(CRISPR_from_prokka_for_all_samples, index=df.index), count_bins

def gather_cas_info_from_prokka(df,annodir,dataset,S3_alternative_dataset_name): # WORKNIG! ma, TODO: da fare TODO hannno megahit nel nome
# WARNING c'è spesso una virgola nel prodigal annotation, dobbiamo liberarcene.
    start_time=time.time()
    print("adding prokka cas entries...")
    os.chdir(annodir+"/justcasanno/")  #questa ancora non è stat divisa in sottocartelle, ma ci mette l'eternita lol
    cas_from_prokka_for_all_samples=[]
    count_bins=0
    for i, row in df.iterrows():
        count_bins+=1
        sample=row["Genome Name"]
        if dataset.startswith("ZeeviD"): #nomenclature exy34yception
            sample=filename_discrepancies.change_to_megahit(sample)
            try:
                sample=sample.replace("ZeeviD_2015", "ZeeviD_2015_B")
                prokkafile_ofsample=open(sample+".cas.gff")
               # print(sample, "uweuwe")
            except:
                sample=sample.replace("ZeeviD_2015_B", "ZeeviD_2015_A")
                prokkafile_ofsample=open(sample+".cas.gff")
               # print(sample, "uweuwe")
        else:
            if dataset != S3_alternative_dataset_name:
                if not (dataset == "ZeeviD_2015_A" or dataset == "ZeeviD_2015_B"):
                    sample=sample.replace(S3_alternative_dataset_name, dataset)
                else:
                    sample, letter=filename_discrepancies.is_this_sample_in_ZeeviD(dataset,sample)
                    if dataset[-1] != letter:
                        continue
            if filename_discrepancies.analysis_has_megahit("prokka"):
                if filename_discrepancies.dataset_has_megahit(dataset, sample):
                    sample_megahit=filename_discrepancies.change_to_megahit(sample)
                    prokkafile_ofsample=open(sample_megahit+".cas.gff")
                else:
#                    print("SONO QUA VA BENE!", os.getcwd(), sample)
                    prokkafile_ofsample=open(sample+".cas.gff")
            else:
                print("This shouldn't even happen")
                break

        cas_from_prokka=""
        for line in prokkafile_ofsample.readlines():
            line=line.rstrip("\n")
            for i in range(len(line)):
                if line[i]==",": # replace the annoying commas inside some annotations with spaces
                    temp=line.split(",")
                    line=" ".join(temp)

            cas_from_prokka +=">"+line.rstrip("\n")
        prokkafile_ofsample.close()
        cas_from_prokka_for_all_samples.append(cas_from_prokka)
    print("Done. Elapsed time: ", str(time.time()-start_time))
    return pd.Series(cas_from_prokka_for_all_samples, index=df.index), count_bins

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
    dataset=sys.argv[1]
    print("++++++++++++++++\nadding CRISPR cas info to dataset: ", dataset)
    outdir="/shares/CIBIO-Storage/CM/scratch/tmp_projects/signorini_cas/3tabellazza/"+dataset
    CRISPRdir="/shares/CIBIO-Storage/CM/scratch/tmp_projects/signorini_cas/1crisprsearch/out"
    annodir="/shares/CIBIO-Storage/CM/scratch/tmp_projects/signorini_cas/2casanno"
    S3_alternative_dataset_name=filename_discrepancies.s3(dataset)
    input_table="/shares/CIBIO-Storage/CM/scratch/tmp_projects/signorini_cas/3tabellazza/" + dataset +"/crisprcas_hits_table_"+dataset+".csv"    #S3Segata_split/S3_subs_"+S3_alternative_dataset_name  #TODO cambia con tabellazze splittate. leva  tutti i punti dall uno al 4 a ggiungi un 5 punto
    print("S3 datasetname = ", S3_alternative_dataset_name)
    df=pd.read_csv(input_table, index_col=0)
    print(df.columns)
    print("initial dataset (rows, cols): ", df.shape)
    print("overall S3 index of current dataset's bins: from ",df.index[0],"to: ", df.index[-1])
#commentate il 22/07 per aggiunta del punto 5. direttamente su crisprcas_hits_table_dataset.csv
#1. add CRISPR loci from the pilercr CRISPR

    pilercr_entry, pilercr_bins=gather_crispr_info_from_PilerCR(df,dataset,CRISPRdir,S3_alternative_dataset_name)
    df["pilercr_CRISPR"]=pilercr_entry
    print(pilercr_entry.shape,  "entries added")
#
##2. add CRISPR loci from the minced CRISPR (if coming from a run with default parameters, these are equal to prokka annotation, therefore from_minced should be set = False)
#
#    if from_minced==True:
#        minced_entry,minced_bins=gather_crispr_info_from_minced(df,dataset, CRISPRdir,S3_alternative_dataset_name)
#        df["minced_CRISPR"]=minced_entry
#        print(minced_entry.shape, " entries added")
#
#
#
#
##3. add CRISPR loci from prokka
#    is_ZeeviD_A, is_ZeeviD_B = False, False
#    if dataset == "ZeeviD_2015_A":
#        is_ZeeviD_A = True
#        print("A")
#    if dataset == "ZeeviD_2015_B":
#        is_ZeeviD_B = True
#        print("B")
#    prokkacrispr,prokka_crispr_bins=gather_crispr_info_from_prokka(df,annodir,dataset,S3_alternative_dataset_name, is_ZeeviD_A, is_ZeeviD_B)
#    print(prokkacrispr.shape, " entries added")
#    df["prokka_CRISPR"]=prokkacrispr
#
##4. add cas loci from prokka
#
#    prokkacas, prokka_cas_bins=gather_cas_info_from_prokka(df,annodir,dataset,S3_alternative_dataset_name)
#    print(prokkacas.shape, " entries added")
#    df["prokka_cas"]=prokkacas

##5. add cas loci from uniref
#
#    unirefcas=gather_cas_info_from_uniref(df,dataset,S3_alternative_dataset_name)
#    print(unirefcas.shape, " entries added")
#    df["uniref_cas"]=unirefcas
#

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
