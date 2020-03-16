import sys
"""USAGE fast_genome_downloader.py fastafile"""
filename=sys.argv[1]
scpcall=""
f=open(filename)
for line in f.readlines():
    if line.startswith(">"):
        line=line.strip().lstrip(">").split("__")
        SGB=line[0]
        dataset=line[2]
        genome=line[2]+"__"+line[3]+"__"+line[4].split(" ")[0]
        mincedcrisprfile="/shares/CIBIO-Storage/CM/scratch/tmp_projects/signorini_cas/1crisprsearch/out/"+dataset+"/"+genome+".fa.minced.out"
        scpcall+="scp lorenzo.signorini@cm3.cibio.unitn.it:"+mincedcrisprfile+" .\n"
f.close()
f=open("scpcommand","w")
f.write(scpcall)
f.close()

