# -*- coding: utf-8 -*-


# Made by L-F-S
# At the University Of Trento, Italy
#
# Define the Cas Locus as a Class, for easier output handling. test
# in vista di un miner veloce che lo usi per spostarsi agilmente tra
# le funzioni


import sys

sys.path.insert(0, '/home/lorenzo.signorini/cas_mining/utils/')
import filename_discrepancies
from Bio.Seq import Seq


outdir="/shares/CIBIO-Storage/CM/news/users/lorenzo.signorini/" #TODO add nstepdir
datadir="/shares/CIBIO-Storage/CM/news/users/lorenzo.signorini/" #TODO add nstepdir

class locus:
    """
        Basic usage:
        import locus

        mylocus=locus.locus(seqid,feature)
        mylocus.fill_from_dataset(cas_dataset)
        tmp=locus.CRISPRarray(feature=feature, contigname=caslocus.contigname,genomename=caslocus.genomename,datasetname=caslocus.datasetname
        tmp.get_CRISPR_array()
             if crarraystrand==-1:
                 tmp.rev_comp()
        caslocus.CRISPRarray=tmp
    """

    def __init__(self, locusname=None, feature=None,seq=None,\
                 CRISPRarray=None,tracrRNA=None, locus_proteins=None,\
                 metadata=None,positions=None, directions=None,\
                 contigname=None, genomename=None, datasetname=None, orig_genomename=None, orig_datasetname=None, estimated_taxonomy=None):
        self.locusname=locusname  # = seqid
        self.feature=feature
        self.seq=seq
        self.CRISPRarray=CRISPRarray
        self.tracrRNA=tracrRNA
        self.locus_proteins=locus_proteins # list of Biopython proteins objects?
        self.metadata=metadata
        self.positions=positions # dictionary of proteins, things, and their start and stop
        self.directions=directions  #dictionary of proteins and things and their direction
        self.contigname=contigname
        self.genomename=genomename
        self.datasetname=datasetname
        self.orig_genomename=orig_genomename
        self.orig_datasetname=orig_datasetname
        self.estimated_taxonomy=estimated_taxonomy


    def fill_from_dataset(self, cas_dataset):
        self.SGB=str(cas_dataset[cas_dataset["Seq ID"]==self.locusname]["SGB ID"].iloc[0])
        self.contigname=str(cas_dataset[cas_dataset["Seq ID"]==self.locusname]["Contig"].iloc[0])
        self.genomename=str(cas_dataset[cas_dataset["Seq ID"]==self.locusname]["Genome Name"].iloc[0])
        self.datasetname=str(cas_dataset[cas_dataset["Seq ID"]==self.locusname]["Study"].iloc[0])
        self.orig_genomename, self.orig_datasetname=filename_discrepancies.get_originalsamplename_froms3name_of_genome(self.genomename,self.datasetname)
        self.seq=cas_dataset[cas_dataset["Seq ID"]==self.locusname]["Seq"].iloc[0]
        self.estimated_taxonomy=cas_dataset[cas_dataset["Seq ID"]==self.locusname]["Estimated taxonomy"].iloc[0]

    def fetch_positions(self, cas_dataset):
        if self.contigname:
            contigname=self.contigname
        else:
            raise ValueError("Please, specify contig name")
        if self.genomename:
            genomename=self.genomename
        else:
            raise ValueError("Please, specify genome name")
        if self.datasetname:
            dataset=self.datasetname
        else:
            raise ValueError("Please, specify dataset name")
        if not self.CRISPRarray.spacers:
            raise ValueError("No spacers found!")
        if not self.CRISPRarray.repeats:
            raise ValueError("No repeats found!")
        if not self.CRISPRarray.repstartpos:
            raise Exception("TODO i feel this is redundant but who knows!")

        cas_position={self.feature:[int(n) for n in cas_dataset[cas_dataset["Seq ID"]==self.locusname]["Pos"].iloc[0].split()],\
                                    "Cas1":[int(n) for n in eval(cas_dataset[cas_dataset["Seq ID"]==self.locusname]["prokka_cas1"].iloc[0])[0][1:]],\
                                    "Cas2":[int(n) for n in eval(cas_dataset[cas_dataset["Seq ID"]==self.locusname]["prokka_cas2"].iloc[0])[0][1:]],\
                                    "CRISPR":[int(n) for n in eval(cas_dataset[cas_dataset["Seq ID"]==self.locusname]["minced_CRISPR"].iloc[0])[0][1:]]
                                   }
        self.positions=cas_position

    def print_locus(self):
        """TODO PRINTS everything tipo output"""




class CRISPRarray:
    def __init__(self, locus=None, feature=None, spacers=None, repeats=None, contigname=None,\
                 genomename=None, datasetname=None,repstartpos=None, path=None):
        self.locus=locus
        self.feature=feature
        self.spacers=spacers
        self.repeats=repeats
        self.contigname=contigname
        self.genomename=genomename
        self.datasetname=datasetname
        self.repstartpos=repstartpos  #TODO  da levar esecodnoeme infuturo
        self.path=path

    def get_CRISPR_array(self, v=True, algorithm="minced"):
        """Retrieve Crispr array position, spacers, and repeats
        algorithm=[\"minced\",\"pilercr\"]
        output: 1) list of spacers, 2) list of repeats, 3) list of start  position of
        every repeat (will be used for BLAST output filtering)"""
        if self.contigname:
            contigname=self.contigname
        else:
            raise ValueError("Please, specify contig name")
        if self.genomename:
            genomename=self.genomename
        else:
            raise ValueError("Please, specify genome name")
        if self.datasetname:
            dataset=self.datasetname
        else:
            raise ValueError("Please, specify dataset name")

        spacers=[]
        repeats=[]
        repeat_start_pos=[]
        if algorithm=="minced":
            orig_genomename, orig_dataset=filename_discrepancies.get_originalsamplename_froms3name_of_genome(genomename,dataset)
            mincedCRISPRfilename=datadir+"1crisprsearch/out/"+dataset+"/"+genomename+".fa.minced.out"
            self.path=mincedCRISPRfilename
            f=open(mincedCRISPRfilename, "r")
            is_the_right_contig=False
            for line in f.readlines():
                if is_the_right_contig:
#                     if v:
#                         print(line)
                    if line[0]== "1" or line[0]== "2" or line[0]== "3" or \
                            line[0]== "4" or line[0]== "5" or line[0]== "6" or\
                            line[0]== "7" or line[0]== "8" or line[0]== "9" or\
                            line[0]== "0":
                        repeats.append(line.split("\t")[2])
                        repeat_start_pos.append(line.split("\t")[0])
                        if not line.split("\t")[3]=="\n":
                            spacers.append(line.split("\t")[3])
                    if line.startswith("Repeats"):
                        break
                if line.startswith("Sequence \'"): #identifiy line of contigname
                    if line.startswith("Sequence \'"+contigname): #check which contig we're on (if there is more than one)
                        is_the_right_contig=True
                    else:
                        is_the_right_contig=False
            f.close()
            self.spacers=spacers
            self.repeats=repeats
            self.repstartpos=repeat_start_pos #TODO da levare secondo me in futuro
            return

        elif algorithm.lower()=="pilercr":
            print("TODO")
        else:
            raise Exception("Allowed algorithms: minced, pilercr")

    def rev_comp(self):
        self.spacers = [str(Seq(spacer).reverse_complement()) for spacer in self.spacers]
        self.repeats = [str(Seq(repeat).reverse_complement()) for repeat in self.repeats]



