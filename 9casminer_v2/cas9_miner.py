#  <h1><center>Data mining for a working Cas9 variant.</center></h1>
# First created Thu Jul 11 09:51:15 CEST 2019
# This is a slightly older version, refer to the Cpf1 pipleine for a more updated version.
# Made by L-F-S
# At the University Of Trento, Italy
import argparse
from os import remove
import sys
sys.path.insert(0, '/home/lorenzo.signorini/cas_mining/9casminer_v2/IO/')
import cas_miner_IO as cmIO

parser=argparse.ArgumentParser(\
                               description="++++++++++++    Extract effector Cas loci from metagenomic data.    ++++++++++++ 01/12/2019: Only extracts Effector Cas sequence :(",\
                               formatter_class=argparse.MetavarTypeHelpFormatter)
#TODO proably dentro il type del parser ci puoi mettere le funzioni tue
#todo add pltting options, es le lengths distrib
parser.add_argument("-v", action="store_true", help="verbose output")
parser.add_argument("-feature", type=str, help="feature name (default: Cas9)", choices=["Cas9","Cpf1"], default="Cas9")
parser.add_argument("-active", type=str, help="Active locus. An active locus is defined as having all proteins and sequences for its feature and type. (default: y)", choices=["y","n"], default="y")
parser.add_argument("-species", type=str, help="species or taxonomiical level (accepted input: [<Genus species>, <any_taxa_name>]) CASE INSENSITIVE")
parser.add_argument("-uSGB", type=str, help="look for unknown (u), known (k) or both (b) SGBs (default: b)", choices=["u","k","b"], default="b")
parser.add_argument("-SGB", type=int, help="look for a specific SGB(s). Space separated")
parser.add_argument("-genome", type=str, help="look for a specific genome(s). Space separated")
parser.add_argument("-length", type=int, nargs="*", help="filter by sequence length, or a range of lengts. Accepted input: <n>, <nmin> <nmax>")
parser.add_argument("-ki", action="store_true", help="FLAG: do not cluster identical sequences, but keep them separated")
parser.add_argument("-o", type=str, help="Set output path. Default: /shares/CIBIO-Storage/CM/scratch/tmp_projects/signorini_cas/5caslocitable/out/", default="/shares/CIBIO-Storage/CM/scratch/tmp_projects/signorini_cas/5caslocitable/out/")
args = parser.parse_args()

def vprint(string,argument=None):
    if argument:
        if args.v:
            print(string)

vprint("Importing modules..",True)
import os
import numpy as np
import pandas as pd
#import matplotlib.pyplot as plt
sys.path.insert(0, '/home/lorenzo.signorini/cas_mining/utils/')
outpath=args.o
datadir="/shares/CIBIO-Storage/CM/scratch/tmp_projects/signorini_cas/"

# # 1. Data:
# ## 1.1 load information about SGBs (species)
vprint("loading data..")
SGB_table=pd.read_csv("/shares/CIBIO-Storage/CM/scratch/tmp_projects/signorini_cas/S4Segata.csv", index_col=0)
# ## 1.2 load information about Cas9s
feature=args.feature
cas_dataset=pd.read_csv("/shares/CIBIO-Storage/CM/scratch/tmp_projects/signorini_cas/5caslocitable/known_"+feature+"_variants_table.csv", index_col=0)

# # 2. Data filtering: Extract a list of the shortest working cas9 from most abundant and most unknown genomes.

# ## 2.1 Filter by: active locus

#Drop Nan rows
if args.active=="y":
    vprint("selecting active loci..", True)
    cas_dataset=cas_dataset.dropna(how="any")
# Subset by species
vprint("Selecting species (or taxonomical level)...", args.species)
cas_dataset=cmIO.subset_by_species(cas_dataset, args.species)
# ## 2.2 Filter by: working active locus
# A working locus is defined by being an active locus with a Cas9 of a length falling within the peak of the distribution
#input("Would you like to plot length distributions?")
#if False:
#    activecas9counts=activecas9s["Seq"].str.count("")
#    plt.figure(figsize=(20*0.8,20*0.8))
#    plt.hist(activecas9counts ,bins=80)
#    plt.xticks(np.arange(0, 1700, step=50))
#    plt.plot()

# The length intreval where there seems to be a working protein is around (1050,1170) , (1330,1450) and maybe (1500,1580). We shall define these intrevals as working lengths.
# filter by SGB:
vprint("Filtering by SGB...", args.SGB)
cas_dataset=cmIO.subset_by_SGB(cas_dataset,args.SGB)
# ## 2.3-2.4 Filter by: sequence length of 949-1099 amino acid, unknown SGB
vprint("Filtering by length...", args.length)
cas_dataset=cmIO.subset_by_length(cas_dataset, args.length)
# filter by unknown SGB An unknown SGB is an SGB with no reference genome in literature.
vprint("Filtering by uSGB...", args.uSGB)
cas_dataset=cmIO.subset_by_unknownSGB(cas_dataset,args.uSGB)#[active_working_cas9s.uSGB=="Yes"]
# filter by genome
vprint("Filtering by genome name...", args.genome)
cas_dataset=cmIO.subset_by_genome(cas_dataset, args.genome)
print("\nDONE!\n\n[TODO:it's actually the total number of genomes,there might be  more than one cas per genome]\n\t--->    Total number of "+args.feature+":\t", len(cas_dataset.index),"  <---",\
      "\n\n- Parameters:\n", args)

# ## show  most abundant SGB

print("Let us check which SGBs are the most abundant among the Cas loci we have filtered so far...")


SGB_abundance_in_dataset=cas_dataset.groupby(["SGB ID"]).count().sort_values(by="Seq ID", ascending=False)["Seq ID"]
SGB_abundance_in_dataset=SGB_abundance_in_dataset.rename("# Genomes") #with act work cas9
SGB_rel_ab_in_dataset=pd.DataFrame(SGB_table[SGB_table["SGB ID"].isin(SGB_abundance_in_dataset.index)][["SGB ID","# Reconstructed genomes"]])
SGB_abundance=pd.DataFrame(SGB_abundance_in_dataset).merge(SGB_rel_ab_in_dataset, left_on='SGB ID', right_on='SGB ID')
SGB_abundance["Genomes relative abundance"]=SGB_abundance["# Genomes"]/SGB_abundance["# Reconstructed genomes"]
print(SGB_abundance)



# # 3 Data Clustering
# ## 3.1 Cluster together identical sequences in the same SGB
# for the 12 most abundant unknown SGBs with a working active Cas9s

# 3. Sequences CLUSTERING
# cluster identical sequences
#todo if args.ki:
print("Extracting "+ feature + "sequences and creating fasta files...")
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO
from Bio.Seq import Seq
from itertools import groupby

#define names for file name:
speciesstr, SGBstr, lengthstr, activestr, genomestr, = "","", "", "", ""
kistr="_identical"
uSGBstr="_"+args.uSGB+"SGB"
if args.species:
    species="_"+args.species
if args.SGB:
    SGBstr="_"+args.SGB+"SGB"
if args.length:
    lengthstr="_".join([str(n) for n in args.length])
if args.active=="y":
    activestr="_active"
if args.genome:
    genomestr="_args.genome"
if args.ki:
    genomestr=""

def seq_getter(s): return str(s.seq)

all_identical_fasta_entries=""

for SGB in np.unique(cas_dataset["SGB ID"]):
    if args.v:
        print("Extracting sequences for SGB: ", SGB)
    current_SGB=cas_dataset[cas_dataset["SGB ID"]==SGB]
    alignments=[]
    for index, row in current_SGB.iterrows():
        tempseq=SeqRecord(Seq(row.Seq), id=str(row["SGB ID"])+"__"+row["Seq ID"]+"__"+row["Genome Name"], description=row["Genome Name"]+"__"+row["Seq ID"]+"__"+str(row["SGB ID"]))
        alignments.append(tempseq)

    #write fasta for every SGB
    filename=outpath+args.feature+"/temp.faa"
    SeqIO.write(alignments, filename, "fasta")
    #read from that fasta and group identical sequences somehow
    if args.v:
        print("Merging Identical sequences for SGB: ", str(SGB))
    records = list(SeqIO.parse(filename,'fasta'))
    records.sort(key=seq_getter)
    n=0
    for seq,equal in groupby(records, seq_getter):
        ids = ';'.join(s.id for s in equal)
        N=ids.count(';')+1

        line=">"+ids+";length"+str(len(seq))+";#sequences"+str(N)+"\n"+seq+"\n"
        all_identical_fasta_entries+=line
    remove(filename)

filename2= outpath+args.feature+"/"+args.feature+"_sequences_"+lengthstr+uSGBstr+SGBstr+activestr+speciesstr+genomestr+kistr+".faa"
f=open(filename2,"w")
f.write(all_identical_fasta_entries)
f.close()
sequences=list(SeqIO.parse(filename2,'fasta'))
print("There are", len(sequences),"unique sequences")
#
for record in  sequences:
    print(record.id)
    print(record.seq)
## ## 3.3 Sequences clustering
## Cluster together sequences of up to 90% similarity, and extract one representative sequence for each cluster.
## Different clustering algorithms are available,the one that works the easiest was fast uclust https://drive5.com/usearch/manual/uclust_algo.html. Uclust performs centroid-based clustering and returns the centroid of each cluster as representative sequence.
##
## ```
## usearch -cluster_fast 949_1099_aa_Cas9_sequences_from_active_loci_from_12_most_abundant_most_unknown_species_identical.faa -id 0.90 -centroids 9aCsfalf1musi_centroid90.faa -uc clusters.uc
##
## ```
#
## uclust output parsing
#
#
#for identity_score in ["90", "93","95","97", "99"]:
#    without_9710, without_8769, without_15286, without_14454, without_9280, without_sketcy_sgbs=[], [],[],[],[],[]
#    n_cluster_per_SGB=[]
#    for record in SeqIO.parse(outpath+identity_score+"_clusters/9aCsfalf1musi_centroid"+identity_score+".faa",'fasta'):
#        if not (record.id.startswith("9710") or record.id.startswith("8769") or record.id.startswith("15286") or record.id.startswith("14454") or record.id.startswith("9280")):
#            without_sketcy_sgbs.append(SeqRecord(record.seq, id=record.id, description=record.description))
#
#        SGBsnames=[oneid[:6].rstrip("__") for oneid in record.id.split(";")]
#        SGBsnames=np.unique(SGBsnames)
#        SGBsnames=",".join(SGBsnames)
#        n_cluster_per_SGB.append(SGBsnames)
#    print("There are", len(n_cluster_per_SGB), "sequence clusters","("+str(len(without_sketcy_sgbs))+" without SGBs 9710,8769, 15286,14454, 9280) at "+identity_score+" identity.")
#    print("SGB(s)\t# Sequence clusters\n---------------------------")
#    print(pd.Series(n_cluster_per_SGB).value_counts(),"\n")
#
#
## ## 3.4 (16/09/2019) MSA withing clusters:
#
#
#identity_score='97'
#clusters="9aCsfalf1musi_centroid"+identity_score+".faa"
#all_sequences_file='949_1099_aa_Cas9_sequences_from_active_loci_from_12_most_abundant_most_unknown_species_identical.faa'
#tot=0
#n=1
#for record in SeqIO.parse(outpath+identity_score+"_clusters/9aCsfalf1musi_centroid"+identity_score+".faa",'fasta'):
#   # collect all sequences from that cluster in this file
#    print("recovering sequences from cluster n", n)
#    alignname="9aCsfalf1musi_clusters"+identity_score+"_all_sequences_from_cluster_"+str(n)
#
#    seqnames=[oneid for oneid in record.id.split(";")]
#    seqnames=seqnames[:-2]
#
#    alignments=[]
#    print(seqnames)
#    for name in seqnames:
#        tot+=1
#        SGB=name[:6].rstrip("__")
#  #      print(SGB)
#        with open(outpath+"SGB"+SGB+".faa", "r") as handle:
#            for record in SeqIO.parse(handle, "fasta"):
#                if name in record.id:
#             #       print("WWWWWWWEEEEEEEEEEEEEEE", name)
#               #     print(record.id)
#                    tempseq=SeqRecord(record.seq, id=record.id)
#                    alignments.append(tempseq)
#    print(len(seqnames),len(alignments))
#
# #   SeqIO.write(alignments, outpath+"97_clusters/"+alignname+".faa", "fasta")
#
#    cline= ClustalwCommandline("clustalw", infile=outpath+"97_clusters/"+alignname+".faa", outfile=outpath+"97_clusters/"+alignname+".aln")
##    print(cline)
#
#    n+=1
#
#print("total amount of sequences= ", tot)
#
#
## # 4 MSalignment with known sequences
##
## perform pairwise alignment of cluster centroids with known workingcas9s. And visualize them nicely.
#
## In[2]:
#
#
## build a fasta file of the things to be aligned against
#
#
## In[71]:
#
#
## merge the fastas together
#identity_score='90'
#clusters="9aCsfalf1musi_centroid"+identity_score+".faa"
#ref_fasta="uniprot_working_Cas9s.fasta"
#alignname="MSA_clusters"+identity_score+"_and_ref"
#
#
## In[74]:
#
#
#alignments=[]
#for row in list(SeqIO.parse(outpath+identity_score+"_clusters/"+clusters,'fasta'))+list(SeqIO.parse(outpath+ref_fasta, 'fasta')):
#    tempseq=SeqRecord(row.seq, id=row.id)
#    alignments.append(tempseq)
#
#SeqIO.write(alignments, outpath+alignname+".faa", "fasta")
#
#
## In[75]:
#
#
#alignname="9aCsfalf1musi_centroid95.faa" #TODO temp
#cline= ClustalwCommandline("clustalw", infile=outpath+alignname+".faa", outfile=outpath+alignname+".aln")
#print(cline)
#os.system(str(cline))
#
#
## In[7]:
#
#
##from Bio import AlignIO
##identity_score='95'
##for record in SeqIO.parse(outpath+"9aCsfalf1musi_centroid"+identity_score+".faa", "fasta"):
# #   for referencecas9record in SeqIO.parse(outpath+"uniprot_working_Cas9s.fasta", "fasta"):
#  #      alignments = pairwise2.align.globalxx(record.seq, referencecas9record.seq)
#        #print(pairwise2.format_alignment(*alignments[0]))
#      #  print("##########################################################################")
#   #     AlignIO.write(alignments, outpath+"test_aln","fasta"),
#
#
#
## # Build maximum likelihood tree out of the multiple sequence alignments
#
## In[1]:
#
#
## go on your laptop open MEGA and build it
#
#
## # Build a score matrix
#
#
#
#from Bio import pairwise2
#import  Bio.SubsMat.MatrixInfo as mx
#
#
## In[7]:
#
#
#help(mx)
#
#
## In[4]:
#
#
##scorematrix=pd.DataFrame()
#scrmtrx=[]
#cols=[]
#for ref_cas9 in SeqIO.parse(outpath+"uniprot_working_Cas9s.fasta", "fasta"):
#    ref_cas9_len=len(ref_cas9.seq)
#    print(">>>>>>>>>>>>>>>>>>>>>>>>REF: "+ref_cas9.id+" ",ref_cas9_len)
#    cols.append(ref_cas9.id)
#    scores=[]
#    index=[] #slight redundancy, this gets called multiple times, ma chissene
#    for mg_cas9 in SeqIO.parse(outpath+"/97_clusters/9aCsfalf1musi_centroid97.faa", "fasta"):
#
#        #aligner = Align.PairwiseAligner()
#        #aligner.substitution_matrix = mx.gonnet
#        #aligner.open_gap_score = -10
#        #aligner.extend_gap_score = -0.1
#        #alignments = aligner.align(mg_cas9.seq, ref_cas9.seq)
#        alignments = pairwise2.align.globalds(mg_cas9.seq, ref_cas9.seq, mx.gonnet, -10, -0.1)
#
#        for ind, (align1, align2, score, begin, end) in enumerate(alignments):
#            filename = outpath+"/97_clusters/pairwise_alignments_of_centroids_vs_refs/PA_" + ref_cas9.id + "_vs_" + mg_cas9.id[:10] + ".aln"
#            with open(filename, "w") as handle:
#                handle.write(">%s\n%s\n>%s\n%s\n" % (mg_cas9.id, align1, ref_cas9.id, align2))
#        print("Done")
#
#        #AlignIO.write(alignments ,outpath+"/97_clusters/pairwise_alignments_of_centroids_vs_refs/" + ref_cas9.id + "_vs_" + mg_cas9.id + ".aln", "clustal")
#
#        score=alignments[0][2]
#        scores.append(score/ref_cas9_len)
# #       scorematrix.at[mg_cas9.id, ref_cas9.id]=score/ref_cas9_len
##        print("SCORE= ", score/ref_cas9_len)
#        index.append(mg_cas9.id)
#  #  scorematrix[ref_cas9.id]=pd.Series(scores)
#    scrmtrx.append(scores)
##scorematrix.index=index
#
#
## In[20]:
#
#
#outpath+"/97_clusters/pairwise_alignments_of_centroids_vs_refs/PA_" + ref_cas9.id + "_vs_" + mg_cas9.id + ".aln"
#
#
## In[15]:
#
#
#help(AlignIO.ClustalIO)
#
#
## In[4]:
#
#
#scrmtrx=np.array(scrmtrx)
#
#
## In[79]:
#
#
##plot a nice heatmap
#import matplotlib.ticker as mtick
#
#
## In[80]:
#
#
#def heatmap(data, row_labels, col_labels, fsz, ax=None,
#            cbar_kw={}, cbarlabel="", **kwargs):
#    """
#    Create a heatmap from a numpy array and two lists of labels.
#
#    Parameters
#    ----------
#    data
#        A 2D numpy array of shape (N, M).
#    row_labels
#        A list or array of length N with the labels for the rows.
#    col_labels
#        A list or array of length M with the labels for the columns.
#    fsz
#        An int to determine font size
#    ax
#        A `matplotlib.axes.Axes` instance to which the heatmap is plotted.  If
#        not provided, use current axes or create a new one.  Optional.
#    cbar_kw
#        A dictionary with arguments to `matplotlib.Figure.colorbar`.  Optional.
#    cbarlabel
#        The label for the colorbar.  Optional.
#    **kwargs
#        All other arguments are forwarded to `imshow`.
#    """
#
#    if not ax:
#        ax = plt.gca()
#
#    # Plot the heatmap
#    im = ax.imshow(data, **kwargs)
#
#    # Create colorbar
#    cbar = ax.figure.colorbar(im, ax=ax, **cbar_kw)
#    cbar.ax.set_ylabel(cbarlabel, rotation=-90, va="bottom", fontsize=fsz)
#
#    # We want to show all ticks...
#    ax.set_xticks(np.arange(data.shape[1]))
#    ax.set_yticks(np.arange(data.shape[0]))
#    # ... and label them with the respective list entries.
#    ax.set_xticklabels(col_labels, fontsize=fsz)
#    ax.set_yticklabels(row_labels, fontsize=fsz)
#
#    # Let the horizontal axes labeling appear on top.
#    ax.tick_params(top=True, bottom=False,
#                   labeltop=True, labelbottom=False)
#
#    # Rotate the tick labels and set their alignment.
#    plt.setp(ax.get_xticklabels(), rotation=-30, ha="right",
#             rotation_mode="anchor")
#
#    # Turn spines off and create white grid.
#    for edge, spine in ax.spines.items():
#        spine.set_visible(False)
#
#    ax.set_xticks(np.arange(data.shape[1]+1)-.5, minor=True)
#    ax.set_yticks(np.arange(data.shape[0]+1)-.5, minor=True)
#    ax.grid(which="minor", color="w", linestyle='-', linewidth=3)
#    ax.tick_params(which="minor", bottom=False, left=False)
#
#    return im, cbar
#
#
#def annotate_heatmap(im, fsz, data=None, valfmt="{x:.2f}",
#                     textcolors=["black", "white"],
#                     threshold=None, **textkw):
#    """
#    A function to annotate a heatmap.
#
#    Parameters
#    ----------
#    im
#        The AxesImage to be labeled.
#    fsz
#        An int to determine font size
#    data
#        Data used to annotate.  If None, the image's data is used.  Optional.
#    valfmt
#        The format of the annotations inside the heatmap.  This should either
#        use the string format method, e.g. "$ {x:.2f}", or be a
#        `matplotlib.ticker.Formatter`.  Optional.
#    textcolors
#        A list or array of two color specifications.  The first is used for
#        values below a threshold, the second for those above.  Optional.
#    threshold
#        Value in data units according to which the colors from textcolors are
#        applied.  If None (the default) uses the middle of the colormap as
#        separation.  Optional.
#    **kwargs
#        All other arguments are forwarded to each call to `text` used to create
#        the text labels.
#    """
#
#    if not isinstance(data, (list, np.ndarray)):
#        data = im.get_array()
#
#    # Normalize the threshold to the images color range.
#    if threshold is not None:
#        threshold = im.norm(threshold)
#    else:
#        threshold = im.norm(data.max())/2.
#
#    # Set default alignment to center, but allow it to be
#    # overwritten by textkw.
#    kw = dict(horizontalalignment="center",
#              verticalalignment="center")
#    kw.update(textkw)
#
#    # Get the formatter in case a string is supplied
#    if isinstance(valfmt, str):
#        valfmt = mtick.StrMethodFormatter(valfmt)
#
#    # Loop over the data and create a `Text` for each "pixel".
#    # Change the text's color depending on the data.
#    texts = []
#    for i in range(data.shape[0]):
#        for j in range(data.shape[1]):
#            kw.update(color="white")#textcolors[int(im.norm(data[i, j]) > threshold)])
#            text = im.axes.text(j, i, valfmt(data[i, j], None), fontsize=fsz, **kw)
#            texts.append(text)
#
#    return texts
#
#
## In[81]:
#
#
#fig, ax = plt.subplots(figsize=(100,50))#300*0.5,200*0.5))
#fontsize=20
#im, cbar = heatmap(scrmtrx, cols, index, fontsize, ax=ax,
#                   cmap="magma", cbarlabel="normalized score")
#texts = annotate_heatmap(im, fontsize, valfmt="{x:.2f}")
#
#fig.tight_layout()
#plt.show()
#
#
## In[49]:
#
#
#fig.savefig("prova.pdf", format="pdf")
#
#
## In[ ]:
#
#
#
#
