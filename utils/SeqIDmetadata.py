# Mon Dec 2 13:27:35 UTC 2019
# Made by L-F-S
# At the University Of Trento, Italy

from pandas import read_csv
import argparse

parser=argparse.ArgumentParser(formatter_class=argparse.MetavarTypeHelpFormatter, description="Print SGB and taxonomical information for a given genome, or feature ID")
parser.add_argument("-feature", type=str,help="Cas9, Cpf1... default: Cas9", default="Cas9")
parser.add_argument("seqid", type=str,help="Effector protein Seq ID")
args=parser.parse_args()

feature=args.feature
outdir="/shares/CIBIO-Storage/CM/scratch/tmp_projects/signorini_cas/Nstep" #TODO add nstepdir
datadir="/shares/CIBIO-Storage/CM/scratch/tmp_projects/signorini_cas/" #TODO add nstepdir
annodir="/shares/CIBIO-Storage/CM/scratch/tmp_projects/signorini_cas/2casanno/crisprcasanno"
cas_dataset=read_csv("/shares/CIBIO-Storage/CM/scratch/tmp_projects/signorini_cas/5caslocitable/known_"+feature+"prokkacas9_variants_table.csv", index_col=0)

seqid=args.seqid
print("SEQ ID: ", seqid)
print("Contig Name:", cas_dataset[cas_dataset["Seq ID"]==seqid]["Contig"].iloc[0])
print("Genome Name:", cas_dataset[cas_dataset["Seq ID"]==seqid]["Genome Name"].iloc[0])
print("SGB: ", cas_dataset[cas_dataset["Seq ID"]==seqid]["SGB ID"].iloc[0])
print("Sequence length: ",len(cas_dataset[cas_dataset["Seq ID"]==seqid]["Seq"].iloc[0]))
print(cas_dataset[cas_dataset["Seq ID"]==seqid]["minced_CRISPR"].iloc[0])
