# Tue May 5 12:08:02 CEST 2020
# Made by L-F-S
# At the University Of Trento, Italy

import os
import subprocess
import argparse
import sys
import time
#import numpy as np
import pandas as pd
from multiprocessing import Pool
import matplotlib.pyplot as plt
from Bio import SeqIO
from Bio.Seq import Seq
#from Bio.SeqRecord import SeqRecord
# for clustal alignments whch i dont use anymore so TODO remove
#from Bio import AlignIO
#from Bio.Align.Applications import ClustalwCommandline
#from Bio.Align import AlignInfo
# to make logos
import logomaker as lm

sys.path.insert(0, '/home/lorenzo.signorini/cas_mining/utils/')
import filename_discrepancies
import originalpath
sys.path.insert(0, '/home/lorenzo.signorini/cas_mining/9casminer_v2/')
import locus
# run on multiple threads
#import multiprocessing as mp
###############################################################################
#               Functions

def vprint(string,v,l):
    if v:
        print(string)
    if l:
      f=open(outdir+"log","a")
      f.write(string+"\n")
      f.close()
    return

def build_blast_db(target_file_name, target_file_dir):
    dbfile=target_file_dir+target_file_name
    if not target_file_name+".nhr" in os.listdir(target_file_dir):
        vprint("Building a blast database for sample: "+ target_file_name+"\n",v,l)
        blastdbcommand="/home/lorenzo.signorini/ncbi-blast-2.10.0+/bin/makeblastdb -in "+dbfile+" -parse_seqids  -dbtype nucl"
        subprocess.Popen(blastdbcommand, shell=True)
    return

def blaster(target_file_name,target_file_dir,outdir,sequid,genomename,feature):
    """ db_file must be a fasta file with
            blast database"""
    dbfile=target_file_dir+target_file_name
    blastoutfile=outdir+seqid+"_vs_"+target_file_name+".blastout"
    blastn_command = "/home/lorenzo.signorini/ncbi-blast-2.10.0+/bin/blastn -out "+blastoutfile+" -outfmt \"6  qseqid sseqid pident qlen length mismatch gapopen qseq sseq sstart send evalue sstrand\" -query temp_spacer_seq -db "+dbfile+" -evalue 0.001 -word_size 12"
    subprocess.Popen(blastn_command, shell=True)
    return


def parallel_wrap_extract_protospacers_and_flanking_sequences(blastoutfile, contigname):
    # TODO logging will be a caso
    temp_samplename=blastoutfile.lstrip(seqid+"_").lstrip("vs_").rstrip("blastout").rstrip(".")
    # check if the code has already been run not in parallel, to save some
    # time:
    if "flanking_sequences_of_putative_protospacers_"+temp_samplename in os.listdir(outdir):
        #vprint("Sample already parsed, moving on to next one.",v,l)
        return
    # run code:
    protospacers_of_sample=extract_protospacers_of(temp_samplename, blastoutfile,contigname)  # returns a dictionary
    if len(protospacers_of_sample.keys())>0:
        extract_flanking_sequences(protospacers_of_sample,temp_samplename) # extrac
    return


def extract_protospacers_of(samplename, blastoutfile,contigname):
    protospacers_of={} # save all putative protospacers position and sequence and samplename, for every contig
    f=open(outdir+blastoutfile)
    lines=f.readlines()
    if not len(lines)>0:
        f.close()
        return {}
    vprint("-"*80+"\n-> Processing sample: "+samplename+"\n",v,l)
    for line in lines:
        target_cont=line.strip("\n").split()[1]  # contig name of target (putative protospacer)
        if target_cont!=contigname: #filter out putative protospacers coming from same contig
            print(line)
            if not target_cont in protospacers_of.keys():
                protospacers_of[target_cont]=[]
            protospacers_of[target_cont].append((int(line.strip("\n").split()[-4]),int(line.strip("\n").split()[-3]),line.strip("\n").split()[-5]))
    if len(protospacers_of.keys())>0:
        vprint("--------> Found "+str(len(protospacers_of.keys()))+" putative protospacers!",v,l)
#    subprocess.Popen("rm "+blastoutfile,shell=True)  # remove blast output file
    return protospacers_of

def extract_flanking_sequences(protospacers_of,target_sample_file):
    """ putative viral contigs extraction
     extract thje +50, -50 sequence flanking protospacer, for every contig
     (ASSUMING 1 contig = 1 protospacer)
     for the particular locus studied, for every contig of putative
    protospacers. Returns None, saves output to file."""
    vprint("Extracting flanking sequences. There are "+str( len(list(protospacers_of.keys())))+  " contigs with putative protospacers in this metagenome.",v,l)
    target_sample_dataset=target_sample_file.split("__")[0]
    dbfile=datadir+target_sample_dataset+"/"+target_sample_file
    f=open(outdir+"flanking_sequences_of_putative_protospacers_"+target_sample_file,"w")
    f=open(outdir+"flanking_sequences_of_putative_protospacers_"+target_sample_file,"a")
    for target_contig in protospacers_of.keys():
        vprint("\n->Processing target contig: "+target_contig+"\n",v,l)

        vprint("Checking if binned...",v,l)
        target_bin='0'
        old_path="/shares/CIBIO-Storage/CM/scratch/users/e.pasolli/projects/binning/genomes_comp50_cont05/"+target_sample_dataset+"/"
        os.chdir(old_path)
        sample_name_only_for_the_following_task=target_sample_file.split(".")[0]
        for filename in os.listdir():
            if filename.startswith(sample_name_only_for_the_following_task):
                for record in SeqIO.parse(filename, "fasta"):
                    if record.id==target_contig:
                        target_bin=filename.rstrip(".fa")
                        vprint("Target contig is in bin: "+target_bin,v,l)
        if target_bin=='0':
            vprint("Target contig is UNBINNED (it was not assigned to any genome).",v,l)

        vprint("Searching for flanking sequences of protospacer from contig: "+target_contig+"\n",v,l)
        for record in SeqIO.parse(dbfile, "fasta"):
            if record.id==target_contig:
                for start, end, matched_seq in protospacers_of[target_contig]:
                    vprint("Protospacer: "+str( matched_seq),v,l)
                    sequence=record.seq
                    revcomp=False
                    if start>end:
                        vprint("Match on minus, taking reverse complement",v,l)
                        s=start
                        start=end
                        end=s
                        revcomp=True
                    test_value=sequence[start-1:end].reverse_complement()==matched_seq if revcomp else sequence[start-1:end] == matched_seq
                    if not test_value:
                        vprint("Exception at sample: "+target_sample_file+",contig: "+target_contig,True,l)
                        raise Exception("Extracted sequence differs:\nExtracted sequence reverse complement?"+str(revcomp))
                    vprint("original protospacer BLAST match:"+str(matched_seq),v,l)
                    vprint("start: "+ str(start)+ "  ,  end: "+str(end),v,l)
                    vprint(str(sequence[start-1:end].reverse_complement()) if revcomp else str(sequence[start:end]),v,l)
                    vprint("extracted protospacer:"+str(test_value),v,l)
                    flanking_length=50
                    if start>=(flanking_length+1):
                        upstream_start=start-(flanking_length+1)
                        added_nucleotides=""
                    else:
                        upstream_start=0
                        added_nucleotides="X"*((flanking_length+1)-start)
                    upseq=added_nucleotides+sequence[upstream_start:start-1] #immediately flanking, no overlap

                    last_pos=len(sequence)
                    if last_pos-end>=flanking_length:
                        downstream_end=end+flanking_length
                        added_nucleotides=""
                    else:
                        downstream_end=last_pos
                        added_nucleotides="X"*(flanking_length-(last_pos-end))
                    downseq=sequence[end:downstream_end]+added_nucleotides  #immediately flanking, no overlap
                    if revcomp:
                        temp=downseq
                        downseq=upseq
                        upseq=temp
                    f.write(target_contig+","+str(upseq)+","+str(downseq)+","+target_sample_file+","+target_bin+","+matched_seq+"\n")
                break
    f.close()
    return

def bulid_sequence_logos(dataset_of_flanking_sequences):
    PAMdata=pd.read_csv(outdir+dataset_of_flanking_sequences, header=None)
    PAMdata=PAMdata.dropna() #TODO temporary

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

def wrapper(seqid,feature,datadir,outdir,castabledir,v,l):
    start_time=time.time()
    # define variables by reading from cas loci table
    #####################################################

    cas_dataset=pd.read_csv(castabledir+"known_"+feature+"_variants_table.csv",\
                            index_col=0)
    dataset=cas_dataset[cas_dataset["Seq ID"]==seqid]["Study"].iloc[0]
    genomename=cas_dataset[cas_dataset["Seq ID"]==seqid]["Genome Name"].iloc[0]
    SGB=cas_dataset[cas_dataset["Seq ID"]==seqid]["SGB ID"].iloc[0]
    contigname=cas_dataset[cas_dataset["Seq ID"]==seqid]["Contig"].iloc[0]
    samplename=cas_dataset[cas_dataset["Seq ID"]==seqid]["Sample Name"].iloc[0]
    full_genome_full_path, annotation_folder=originalpath.print_path(genomename)
    # needed for current situation of original data
    old_genomename, old_dataset=filename_discrepancies.get_originalsamplename_froms3name_of_genome(genomename,dataset)
    old_path = "/shares/CIBIO-Storage/CM/scratch/users/e.pasolli/projects/binning/genomes_comp50_cont05/"+old_dataset+"/"
    #
    if not os.path.exists(outdir):
            os.makedirs(outdir)
    # initialize log:
    if l:
        f=open(outdir+"log","w")
        f.close()
    vprint("-"*80,True,l)
    vprint("Sequence ID: "+ seqid,True,l)
    vprint("Contig:      "+ contigname,True,l)
    vprint("Dataset:     "+ dataset,True,l)
    vprint("Genome name: "+ genomename,True,l)
    vprint("Sample name: "+ samplename,True,l)
    vprint("SGB ID:      "+ str(SGB),True,l)
    vprint("Cas9 length: "+str(len(cas_dataset[cas_dataset["Seq ID"]==seqid]["Seq"].iloc[0])),True,l)
    vprint("old epasolli genomename:  "+ old_genomename,True,l)

    # get spacers from CRISPRarray sequence
    ############################################################

    cr=locus.CRISPRarray(feature=feature, contigname=contigname, genomename=genomename, datasetname=dataset)
    cr.get_CRISPR_array()
    spacers=cr.spacers
    vprint("spacers:   "+str(spacers),v,l)

    # build temporary query file for blastn search, with spacers
    ########################################################

    vprint("-"*80+"\n",True,l)
    vprint("building query file of spacers.",True,l)
    blast_folder=outdir
    os.chdir(blast_folder)
    tempfile=open("temp_spacer_seq", "w")
    tempfile.close()
    tempfile=open("temp_spacer_seq", "a")
    for n, spacer in enumerate(spacers):
            tempfile.write(">spacer"+str(n+1)+"|"+contigname+"|"+genomename+"|"+seqid+"\n"+spacer+"\n")
    tempfile.close()
    vprint("time passed: "+ str(time.time()-start_time),v,l)

    # Blast query file against a set of databases
    ############################################################

    #vprint("\nRunning blastn of query file agianst all contigs",True,l)
#   # #TODO add optional user input to set blast params and directiory
    #old_dataset_names=['HeQ_2017', 'BengtssonPalmeJ_2015', 'CM_caritro','LoombaR_2017', 'GeversD_2014', 'QinN_2014', 'SchirmerM_2016', 'RampelliS_2015', 'LiSS_2016', 'HMP_2012', 'QinJ_2012', 'LiJ_2014', 'BritoIL_2016', 'LouisS_2016', 'ZeeviD_2015_B', 'CM_madagascar', 'KosticAD_2015', 'LawrenceA_2015', 'SmitsSA_2017', 'OlmMR_2017', 'IjazUZ_2017', 'VincentC_2016', 'Obregon-TitoAJ_2015', 'CosteaPI_2017', 'KarlssonFH_2013', 'VogtmannE_2016', 'NielsenHB_2014', 'VatanenT_2016', 'YuJ_2015', 'ChngKR_2016', 'ChengpingW_2017', 'ZeeviD_2015_A', 'CM_periimplantitis', 'AsnicarF_2017', 'RaymondF_2016', 'OhJ_2014', 'Castro-NallarE_2015', 'XieH_2016', 'LomanNJ_2013', 'LeChatelierE_2013', 'ZellerG_2014', 'HanniganGD_2017', 'WenC_2017', 'BackhedF_2015', 'CM_cf', 'LiuW_2016', 'LiJ_2017', 'FengQ_2015']
    #for dataset in old_dataset_names:
    #    vprint("Blasting against samples of dataset "+dataset+"..",True,l)
    #    target_dir=datadir+dataset+"/"
    #    for tempdatabase in os.listdir(target_dir):
    #        if tempdatabase.endswith(".fasta"):
    #            target_sample_file=tempdatabase
    #            build_blast_db(target_sample_file, target_dir)
    #            blaster(target_sample_file,target_dir,outdir,seqid,genomename,feature)
    #vprint("Time passed: "+str(time.time()-start_time),True,l)

    # parse BLAST output for putative protospacers &
    # extact flanking regions
    ###################################################################
    #vprint("-"*80,True,l)
    #vprint("Parsing Blast output:\nFinding putative protospacers and flanking regions",True,l)

    # Multi-Process version of for cycle below. Uncomment to execute with TODO
    # 30 cores:

    #alternative input list:
    #inputs_list=[]
    #f=open(datadir+"samples_to_rerun.csv")
    #for line in f.readlines:
    #    inputs_list.append((seqid+"_vs_"+line.strip("\n")+".blastout",contigname))
    #f.close()

    inputs_list=[(blastoutfile,contigname) for blastoutfile in  os.listdir(outdir) if blastoutfile.endswith("blastout")] # create a list of tuples (target_blast_output_file,contigname) to pass to the multiprocess iterator
    pool=Pool(30)
    pool.starmap(parallel_wrap_extract_protospacers_and_flanking_sequences,inputs_list)
    # the rest  will be executed after multiprocesses are completed

    # uncomment to execute with 1 core:
   # for blastoutfile in os.listdir(outdir):  #cycle though all samples
   #     # parse all blast outputs and extract putative protospacers
   #     if blastoutfile.endswith("blastout"):
   #         temp_samplename=blastoutfile.lstrip(seqid+"_").lstrip("vs_").rstrip("blastout").rstrip(".")
   #         protospacers_of_sample=extract_protospacers_of(temp_samplename, blastoutfile,contigname)  # returns a dictionary
   #         if len(protospacers_of_sample.keys())>0:
   #             extract_flanking_sequences(protospacers_of_sample,temp_samplename) # extract flanking sequences from original genomes and saves to a file TODO an easier version would extract from sample file instead of old genomes!!
   ##################################################################

    # Merge all files in one dataset
    vprint("Merging all flanking sequences in one dataset...",True,l)
    os.chdir(outdir)
    dataset_of_flanking_sequences="dataset_flanking_sequences_of_putative_protospacers"
    subprocess.Popen("cat flanking_sequences* > dataset_flanking_sequences_of_putative_protospacers", shell=True)
    #subprocess.Popen("cat flanking_sequences* > dataset_flanking_sequences_of_putative_protospacers && rm flanking_sequences*", shell=True)
   # vprint("Time passed: "+str(time.time()-start_time),True,l)

    # Build sequence logos
    ##################################################################
    vprint("Building upstream and downstream sequence logos...",True,l)
    #TODO add input to modify length of seq logos
    bulid_sequence_logos(dataset_of_flanking_sequences)
    return







##############################################################################
if __name__=="__main__":
    main_descr="Wellcome to PAM_finder! The in-silico PAM sequence discovery essay for metagenomic data! Insert a locus id (effector id). A log file at the end of the process will be in the output folder"
    parser=argparse.ArgumentParser(description="+"*5+"\t\t"+main_descr+"\t\t"+"+"*5 )
    parser.add_argument("ID", type=str, help="sequence ID")
    parser.add_argument("-v", action="store_true", help="verbose output")
    parser.add_argument("-l", action="store_false", help="Do not log PAM_finder's outptut.")
    parser.add_argument("-f", type=str, help="effector Cas name (default= Cas9)"\
                        , default="Cas9")
    parser.add_argument("-c", type=str, help="cas_dataset position, default=\
                        /shares/CIBIO-Storage/CM/news/users/lorenzo.signorini/5caslocitable"\
                        ,default="/shares/CIBIO-Storage/CM/news/users/lorenzo.signorini/5caslocitable/")
    parser.add_argument("-o", type=str, help="output directory, default =\
                        /shares/CIBIO-Storage/CM/news/users/lorenzo.signorini\
                        /8pamsearch/out/", default="/shares/CIBIO-Storage/CM/news/users/lorenzo.signorini/8pamsearch/out/")
    parser.add_argument("-w", type=str, help="data directory, where your contigs\
                        BLAST databases are, default \
                        =/shares/CIBIO-Storage/CM/news/users/lorenzo.signorini/"\
                        , default="/shares/CIBIO-Storage/CM/news/users/lorenzo.signorini/8pamsearch/")
    parser.add_argument("-d", action="store_true",  help="Save sequence output")
    args=parser.parse_args()
    seqid =args.ID
    feature=args.f
    outdir=args.o+seqid+"/"
    datadir=args.w
    castabledir=args.c
    v=args.v
    l=args.l
    if not os.path.exists(outdir):
            os.makedirs(outdir)
    wrapper(seqid,feature,datadir,outdir,castabledir,v,l)
