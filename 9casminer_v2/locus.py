

# Made by L-F-S
# At the University Of Trento, Italy
#
# Define the Cas Locus as a Class, for easier output handling. test
# in vista di un miner veloce che lo usi per spostarsi agilmente tra
# le funzioni


import sys

sys.path.insert(0, '/home/lorenzo.signorini/cas_mining/utils/')
import filename_discrepancies


outdir="/shares/CIBIO-Storage/CM/scratch/tmp_projects/signorini_cas/Nstep" #TODO add nstepdir
datadir="/shares/CIBIO-Storage/CM/scratch/tmp_projects/signorini_cas/" #TODO add nstepdir
annodir="/shares/CIBIO-Storage/CM/scratch/tmp_projects/signorini_cas/2casanno/crisprcasanno"

class locus:

    def __init__(self, locusname=None, feature=None,seq=None,\
                 CRISPRarray=None,tracrRNA=None, locus_proteins=None,\
                 metadata=None,positions=None, directions=None,\
                 contigname=None, genomename=None, datasetname=None):
        self.locusname=locusname
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

    def set_SGB(self, cas_dataset):
        """return SGB as a string. cas_dataset is  5cas_loci_table formatted table"""
        self.SGB=str(cas_dataset[cas_dataset["Seq ID"]==self.locusname]["SGB ID"].iloc[0])

    def set_contigname(self, cas_dataset):
        """return SGB as a string. cas_dataset is  5cas_loci_table formatted table"""
        self.contigname=str(cas_dataset[cas_dataset["Seq ID"]==self.locusname]["Contig"].iloc[0])

    def set_genomename(self, cas_dataset):
        """return SGB as a string. cas_dataset is  5cas_loci_table formatted table"""
        self.genomename=str(cas_dataset[cas_dataset["Seq ID"]==self.locusname]["Genome Name"].iloc[0])

    def set_datasetname(self, cas_dataset):
        """return SGB as a string. cas_dataset is  5cas_loci_table formatted table"""
        self.datasetname=str(cas_dataset[cas_dataset["Seq ID"]==self.locusname]["Study"].iloc[0])
    def set_sequence(self, cas_dataset):
        self.seq=cas_dataset[cas_dataset["Seq ID"]==self.locusname]["Seq"].iloc[0]

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

#    def BLAST_antirepeats(self):
#        genomefilename= "/shares/CIBIO-Storage/CM/scratch/tmp_projects/epasolli_darkmatter/allcontigs/ALLreconstructedgenomes/"+str(SGB)+"/"+genomename+".fa"
#        print("+"*80)
#        print(cas_position)
#        print("Retrieving locus")
#
#
#        # access contig:
#
#        for record in SeqIO.parse(genomefilename,"fasta"):  #c'è modo di non fare questo 'ciclo'?
#            if record.id.startswith(contigname):
#                   #test se effettivamente la seq di cas è dove voglio io
#
#                SeqIO.write(record, "/shares/CIBIO-Storage/CM/scratch/tmp_projects/signorini_cas/7tracrRNA/tempblastdb", "fasta")
#                print("CRISPR si trova", cas_position["CRISPR"])
#                print("Cas9 si trova", cas_position[feature])
#                print("Cas2 si trova", cas_position["Cas2"])
#                print("Cas1 si trova", cas_position["Cas1"])
#
#                 #estrai la posizione di cas9 dal contig, e traducila
#                temp_alt_cas9start=98503
#                temp_alt_cas9stop=101700
#                cas9start=cas_position[feature][0]
#                cas9stop=cas_position[feature][1]
#                cas9_nnseq=record.seq[cas9start-1:cas9stop-3] #gff should have a 1-based positional annotation (ma sto andando a occhio finchè non sono uguali, è già un oggetto Bio.Seq
#                my_translated_cas9=cas9_nnseq.transcribe().translate()
#                print(len(cas9_nnseq)/3,len(cas9_aa))
#                if my_translated_cas9==cas9_aa:
#                    print("Cas locus on plus strand")
#                    plus=True
#                if not my_translated_cas9==cas9_aa:
#                  #try reverse complement
#                    cas9_nnseq=record.seq[cas9start+2:cas9stop]
#                    my_translated_cas9=cas9_nnseq.reverse_complement().transcribe().translate()
#                    print("Cas locus on minus strand")
#                print("Is Cas9 translated from the genome equal to the orignal annotation?", my_translated_cas9==cas9_aa)
#                print(my_translated_cas9)
#                print("+"*80)
#                print(cas9_aa)
#                alignments = pairwise2.align.localxx(my_translated_cas9,cas9_aa)
#                    print(pairwise2.format_alignment(*alignments[0]))
#
#                            #opzione 1 blast
#                                  #1. salva np.uniq(repeats) in un fasta file: build temporary query file for blastn search
#
#                                          print("let's now blast the repeat against  the genome")
#                                                  blast_folder="/shares/CIBIO-Storage/CM/scratch/tmp_projects/signorini_cas/7tracrRNA/"
#                                                          os.chdir(blast_folder)
#                                                                  tempfile=open("temp_repeat_seq", "w")
#                                                                          tempfile.close()
#                                                                                  tempfile=open("temp_repeat_seq", "a")
#                                                                                          for n, repeat in enumerate([Seq(sequence) for sequence in np.unique(repeats)]):
#                                                                                                          tempfile.write(">rpt"+str(n+1)+"|"+contigname+"|"+genomename+"|"+seqid+"\n"+str(repeat)+"\n")
#                                                                                                                  tempfile.close()
#
#                                                                                                                          #2. Blast query file against db
#                                                                                                                                  print("Running blastn of query file agianst all contigs")
#
#                                                                                                                                          #3. make db file
#                                                                                                                                                  dbfile="/shares/CIBIO-Storage/CM/scratch/tmp_projects/signorini_cas/7tracrRNA/tempblastdb"
#                                                                                                                                                          os.system("makeblastdb -in "+dbfile+" -parse_seqids  -dbtype nucl")
#
#
#
#                                                                                                                                                                  blastoutfile="/shares/CIBIO-Storage/CM/scratch/tmp_projects/signorini_cas/7tracrRNA/"\
#                                                                                                                                                                              +str(SGB)+"__"+seqid+"__"+genomename+"__"+feature+".trcr.blastout"
#                                                                                                                                                                          blastn_command = "blastn -out "+blastoutfile+" -outfmt \"6  qseqid sseqid pident qlen length mismatch gapopen qseq sseq sstart send evalue sstrand\" -query temp_repeat_seq -db "+dbfile+" -evalue 0.001 -word_size 11 -penalty -2"
#                                                                                                                                                                                  print(blastn_command)
#                                                                                                                                                                                          os.system(blastn_command)
#                                                                                                                                                                                                 # os.system("rm
#                                                                                                                                                                                                 # /shares/CIBIO-Storage/CM/scratch/tmp_projects/signorini_cas/7tracerRNA/temp*")
#



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
                    if v:
                        print(line)
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


