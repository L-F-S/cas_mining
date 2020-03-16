# Wed Mar 11 13:39:44 CET 2020
# Made by L-F-S
# At the University Of Trento, Italy

import sys



outdir="/shares/CIBIO-Storage/CM/scratch/tmp_projects/signorini_cas/Nstep" #TODO add nstepdir
datadir="/shares/CIBIO-Storage/CM/scratch/tmp_projects/signorini_cas/" #TODO add nstepdir
annodir="/shares/CIBIO-Storage/CM/scratch/tmp_projects/signorini_cas/2casanno/crisprcasanno"

"""USAGE fast_genome_downloader.py fastafile"""
filename=sys.argv[1]
scpcall=""
f=open(filename)
for line in f.readlines():
    if line.startswith(">"):
        line=line.strip().lstrip(">").split("__")
        SGB=line[0]
        genome=line[2]+"__"+line[3]+"__"+line[4].split(" ")[0]
        scpcall+="scp lorenzo.signorini@cm3.cibio.unitn.it:/scratchCM/tmp_projects/epasolli_darkmatter/allcontigs/ALLreconstructedgenomes/"+SGB+"/"+genome+".fa .\n"
f.close()
f=open("scpcommand","w")
f.write(scpcall)
f.close()
