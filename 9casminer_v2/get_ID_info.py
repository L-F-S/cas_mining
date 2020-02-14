
# Made by L-F-S

# At the University Of Trento, Italy
# lorenzo.signorini@alumni.unitn.it

"""
For usage information:


python get_ID_info.py -h


"""
#TODO onc eI automatize tracrfinding, we will be able to make this into a
#uuseful tool
import os
import sys
import pandas as pd
import argparse
from Bio import SeqIO

import locus
sys.path.insert(0, '/home/lorenzo.signorini/cas_mining/utils/')
import filename_discrepancies

def get_ID_info(seqid, feature,v, saveout, outdir,tracrRNA, tracrstrand, crarraystrand,repeat):
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

    if v:
        print("\n-> Sequence:\n\n", caslocus.seq)
    tmp=locus.CRISPRarray(feature=feature, contigname=caslocus.contigname,genomename=caslocus.genomename,datasetname=caslocus.datasetname)
    tmp.get_CRISPR_array()
    if crarraystrand==-1:
        tmp.rev_comp()
    #todo da cambiare la riga sopra di me e aggiungere un metodo .print
    #e una line qua: if v: printCIRSPRarray, e levare quel verboso la
    caslocus.CRISPRarray=tmp #Magari slightly redundant to have that locus
    caslocus.fetch_positions(cas_dataset)
    print("-> Positions:\t",caslocus.positions)
    #print("TODO: Finding tracr_RNA")
    caslocus.tracrRNA=tracrRNA
    if args.t:
        print("-> tracrRNA length:" ,len(caslocus.tracrRNA))

    if v:
# fai qui CRISPR
        if args.t:
            print("\n----------------------------------------------------------------------")
            print("         tracrRNA sequence for", feature,"protein id", seqid)
            print("----------------------------------------------------------------------\n")
            print(caslocus.tracrRNA)
        if repeat:
            print("Matching repeat:\n"+repeat)

        CRISPRheader="REPEAT"+" "*(len(caslocus.CRISPRarray.repeats[0])-len("REPEAT"))+" SPACER"+" "*(len(caslocus.CRISPRarray.spacers)-len("SPACERS"))+"\n"
        print("CRISPR array sequence for "+ seqid+" "+feature+".")
        print("CRISPR array on REVERSE strand.")
        print(CRISPRheader)
        CRISPRarraysequence=""
        for i in range(len(caslocus.CRISPRarray.repeats)):
            try:
                caslocus.CRISPRarray.spacers[i]
                line=caslocus.CRISPRarray.repeats[i]+" "+caslocus.CRISPRarray.spacers[i]
                print(line)
                CRISPRarraysequence+=line+"\n"
            except:
                break
    #TODO add the written output
    if saveout:
        print("\n----------------------------------------------------------------------")
        print("Saving output in: \n"+outdir)
        print("----------------------------------------------------------------------\n")
        if not os.path.exists(outdir+seqid):
                os.makedirs(outdir+seqid)
        print("Writing "+feature+" amino acid sequence to "+feature+"_"+seqid+".faa")
        f=open(outdir+seqid+"/"+feature+"_"+seqid+".faa","w")
        f.write(">"+seqid+"\n"+caslocus.seq)
        f.close()
        if args.t:
            print("Writing tracrRNA to tracrRNA_"+seqid+".ffn")
            f=open(outdir+seqid+"/tracrRNA_"+seqid+".ffn","w")
            f.write(">tracrRNA_"+seqid+"\n"+tracrRNA)
            f.close()
        print("Writing  CRISPRarray to CRISPR_"+seqid+".ffn")
        if crarraystrand==-1:
            f=open(outdir+seqid+"/CRISPR_"+seqid+".ffn", "w")
            f.close()
            f=open(outdir+seqid+"/CRISPR_"+seqid+".ffn", "a")
            f.write("@made by L-F-S\n@ University of Trento\n@lorenzo.signorini@alumni.unitn.it\n\n")
            f.write("CRISPR array sequence for "+ seqid+" "+feature+".\nCRISPR on REVERSE strand\n\n\n")
            f.write(CRISPRheader)
            f.write(CRISPRarraysequence)
            os.popen("cp "+tmp.path+" "+outdir+seqid+"/NOT_REAL_CRISPR_but_original_annotation_of_CRISPR_"+seqid+".ffn")
        else:
            os.popen("cp "+tmp.path+" "+outdir+seqid+"/CRISPR_"+seqid+".ffn")



        #V2:
        print("Writing"+feature+" nucleotidic sequence to "+feature+"_seqid"+".ffn")
        #chiaramente appiccicato da un altro file
        genomename, dataset=filename_discrepancies.get_originalsamplename_froms3name_of_genome(caslocus.genomename,caslocus.datasetname)
        s3_genomename=caslocus.genomename
        # qui genomename e datasetname si riferiscono al nome dentro epasolli
        # darkmatter (Zeevid_A, ZeeviD_B, megahit)
        # mamma mia the redundancy of those two lines, vabbe
        contigname=caslocus.contigname
        SGB=caslocus.SGB
        cosa=".ffn"
        prokka_anno_file="/scratchCM/tmp_projects/epasolli_darkmatter/allcontigs/"+dataset+"/metabat/genomes_comp50_cont05/prokka/"+genomename+"/"+genomename+cosa
        for record in SeqIO.parse(prokka_anno_file,"fasta"):
            if record.id.startswith(seqid):
                record.description=feature+" len="+str(len(record.seq))+" genome="+s3_genomename+" SGB="+str(SGB)+" contig="+contigname
                Cas9_fasta_header=">"+seqid+" "+record.description
    #            coso_non_capisco_piu_nulla[cosa][feature]=record.seq
                SeqIO.write(record, outdir+seqid+"/"+feature+"_"+seqid+cosa,"fasta")
    #        if record.id.startswith(Cas2ID):
    #            record.description="Cas2 len="+str(len(record.seq))+" genome="+s3_genomename+" SGB="+str(SGB)+" contig="+contigname
    #            Cas2_fasta_header=">"+Cas2ID+" "+record.description
    #            SeqIO.write(record, outputdir+"Cas2_"+seqid+cosa,"fasta")
    #            coso_non_capisco_piu_nulla[cosa]["Cas2"]=record.seq
    #        if record.id.startswith(Cas1ID):
    #            record.description="Cas1 len="+str(len(record.seq))+" genome="+s3_genomename+" SGB="+str(SGB)+" contig="+contigname
    #            Cas1_fasta_header=">"+Cas1ID+" "+record.description
    #            coso_non_capisco_piu_nulla[cosa]["Cas1"]=record.seq
    #            SeqIO.write(record, outputdir+"Cas1_"+seqid+cosa,"fasta")





if __name__=="__main__":
    main_descr="Print information about a given effector Cas Sequence ID."
    plus="+"*len(main_descr)
    parser=argparse.ArgumentParser(description="+"*5+"\t\t"+main_descr+"\t\t"+"+"*5 )
    parser.add_argument("-v", action="store_true", help="verbose output")
    parser.add_argument("ID", type=str, help="sequence ID")
    parser.add_argument("-f", type=str, help="effector Cas name (default= Cas9)", default="Cas9")
    parser.add_argument("-o", type=str, help="output directory, default =/shares/CIBIO-Storage/CM/scratch/tmp_projects/signorini_cas/9output", default="/shares/CIBIO-Storage/CM/scratch/tmp_projects/signorini_cas/9output/")
    parser.add_argument("-d", action="store_true",  help="Save sequence output")
    parser.add_argument("-t", type=str, help="tracRNA sequence")
    parser.add_argument("-s", type=str, help="tracrRNA strand. Possible values= +, -")
    parser.add_argument("-c", type=int, help="actual CRISPRarray strand (as discovered via previous tracrRNA analysis)\
                        . Input \'+1\' if CRISPRarray strand is the same annotated from minced algorithm (forward strand), \'-1\' if reverse strand.")
    parser.add_argument("-r", type=str, help="repeat sequence matching tracrRNA")
    args=parser.parse_args()
    outdir=args.o
    seqid =args.ID
    feature=args.f
    tracrRNA=args.t
    tracrstrand=args.s
    crarraystrand=args.c
    repeat=args.r

    get_ID_info(seqid, feature,args.v,args.d,outdir,tracrRNA, tracrstrand, crarraystrand,repeat)
