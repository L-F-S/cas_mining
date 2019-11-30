# Fri Nov 29 11:23:26 CET 2019
# Made by L-F-S
# At the University Of Trento, Italy
#
# Define the Cas Locus as a Class, for easier output handling. test
# in vista di un miner veloce che lo usi per spostarsi agilmente tra
# le funzioni


import sys

sys.path.insert(0, '/home/lorenzo.signorini/cas_mining/utils/')
import filename_discrepancies


dataset="SmitsSA_2017" #sys.argv[1] #TODO switch test dataset
outdir="/shares/CIBIO-Storage/CM/scratch/tmp_projects/signorini_cas/Nstep" #TODO add nstepdir
datadir="/shares/CIBIO-Storage/CM/scratch/tmp_projects/signorini_cas/" #TODO add nstepdir
annodir="/shares/CIBIO-Storage/CM/scratch/tmp_projects/signorini_cas/2casanno/crisprcasanno"

class locus:
    def __init__(self, locusname, effector_cas, CRISPRarray, tracrRNA,\
                 locus_proteins, metadata, positions, directions ):
        self.locusname=locusname
        self.effector_cas=effector_cas
        self.CRISPRarray=CRISPRarray
        self.tracrRNA=tracrRNA
        self.locus_proteins=locus_proteins # list of Biopython proteins objects?
        self.metadata=metadata
        self.positions=positions # dictionary of proteins, things, and their start and stop
        self.directions=directions  #dictionary of proteins and things and their direction

    def SGB(self, metadata):
        """return SGB as a string"""
        SGB="TODO"
        return SGB
    def pyhlo(self, metadata):
        """returns S4tab phylogeny"""
        phylo="todo"
        return phylo
    def species(self, metadata):
        """returns genus and species"""
        sp="TODO"
        return sp

    def print_locus(self):
        """TODO PRINTS everything tipo output"""


class CRISPR:
    def __init__(self, name, CRISPRarray):
        self.name=name
        self.CRISPRarray=CRISPRarray

    def get_spacers(CRISPRarray):
        """TODO"""
        spacers="todo"
        return spacers

    def get_repeats(CRISPRarray):
        """TODO"""
        repeats="todo"
        return repeats

    def print_array(CRISPRarray):
        """TODO"""
        return


class tracrRNA:
    def __init__(self,antirepeat,terminator, sequence):
        self.antirepeat=antirepeat
        self.terminator=terminator
        self.sequence=sequence
