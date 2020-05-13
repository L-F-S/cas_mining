
# Made by L-F-S

# At the University Of Trento, Italy
# lorenzo.signorini@alumni.unitn.it

"""
For usage information:


python get_ID_info.py -h


"""
import os
import sys
import pandas as pd
import argparse
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Align.Applications import ClustalwCommandline

import locus
sys.path.insert(0, '/home/lorenzo.signorini/cas_mining/utils/')
import originalpath

def vprint(string):
    if args.v:
        print(string)


def get_ID_info(seqid, feature,v, saveout, outdir,tracrRNA, repeat,wdir):
    os.chdir(wdir)
    # intialize Cas locus object, from cas dataset
    cas_dataset=pd.read_csv(wdir+"5caslocitable/known_"+feature+"_variants_table.csv", index_col=0)
    caslocus=locus.locus(seqid,feature)
    caslocus.fill_from_dataset(cas_dataset)
    # get CRISPRarray and positions
    tmp=locus.CRISPRarray(feature=feature, contigname=caslocus.contigname,genomename=caslocus.genomename,datasetname=caslocus.datasetname)
    tmp.get_CRISPR_array()
    if args.c==-1:
        tmp.rev_comp()
    caslocus.CRISPRarray=tmp
    caslocus.fetch_positions(cas_dataset)
    # Print main info
    info_text=""
    temp_print="-"*80+"\n"+" "*int((80- len("Sequence ID"))/2)+" Sequence ID"+\
          " "*int((80- len("Sequence ID"))/2)+"\n"+\
          " "*int((80- len("Sequence ID"))/2)+ seqid+\
          " "*int((80- len("Sequence ID"))/2)+"\n"+"-"*80+"\n"
    info_text+=temp_print+"\n"
    print(temp_print)
    temp_print="-> Contig Name:\t"+caslocus.contigname
    info_text+=temp_print+"\n"
    print(temp_print)
    temp_print="-> Genome Name:\t"+caslocus.genomename
    info_text+=temp_print+"\n"
    print(temp_print)
    temp_print="-> Dataset Name:\t"+caslocus.datasetname
    info_text+=temp_print+"\n"
    print(temp_print)
    temp_print="-> Epasolli Genome Name: "+ caslocus.orig_genomename
    info_text+=temp_print+"\n"
    print(temp_print)
    temp_print="-> Epasolli Dataset Name: "+ caslocus.orig_datasetname
    info_text+=temp_print+"\n"
    print(temp_print)
    temp_print="-> SGB ID:\t"+caslocus.SGB
    info_text+=temp_print+"\n"
    print(temp_print)
    temp_print="-> Estimated taxonomy:\t"+caslocus.estimated_taxonomy
    info_text+=temp_print+"\n"
    print(temp_print)
    temp_print="-> " +feature+" sequence length:\t"+str(len(caslocus.seq))
    info_text+=temp_print+"\n"
    print(temp_print)
    if v:
       print("\n-> Sequence:\n\n", caslocus.seq)
    temp_print="-> Positions inside contig:\t"+str(caslocus.positions)
    print(temp_print)
    info_text+=temp_print+"\n"
    caslocus.tracrRNA=tracrRNA
    if args.t:
        temp_print="-> tracrRNA length:\t"+str(len(caslocus.tracrRNA))
        info_text+=temp_print+"\n"
        print(temp_print)
        temp_print="-> tracrRNA strand:\t"+args.s
        info_text+=temp_print+"\n"
        print(temp_print)
        temp_print="-> CRISPRarray strand:\t"+str(args.c)
        info_text+=temp_print+"\n"
    CRISPRheader="REPEAT"+" "*(len(caslocus.CRISPRarray.repeats[0])-len("REPEAT"))+" SPACER"+" "*(len(caslocus.CRISPRarray.spacers)-len("SPACERS"))+"\n"
    if args.c==1:
        CRISPRstrandstring="CRISPR array on FORWARD strand"
    elif args.c==-1:
        CRISPRstrandstring="CRISPR array on REVERSE strand."
    else:
        CRISPRstrandstring="Real CRISPR strand unknown, printing annotation in forward direction."
    # verbose output: print CRISPR ARRAY
    vprint("CRISPR array sequence for "+ seqid+" "+feature+".")
    vprint(CRISPRstrandstring)
    vprint(CRISPRheader)
    CRISPRarraysequence=""
    for i in range(len(caslocus.CRISPRarray.repeats)):
        try:
            caslocus.CRISPRarray.spacers[i]
            line=caslocus.CRISPRarray.repeats[i]+" "+caslocus.CRISPRarray.spacers[i]
            vprint(line)
            CRISPRarraysequence+=line+"\n"
        except:
            line=caslocus.CRISPRarray.repeats[i]
            vprint(line)
            CRISPRarraysequence+=line+"\n"
            break
        # verbose output: print tracrRNA sequences
    if args.t:
        vprint("\n----------------------------------------------------------------------")
        vprint("         tracrRNA sequence for "+ feature +" protein id "+ seqid)
        vprint("----------------------------------------------------------------------\n")
        vprint(caslocus.tracrRNA)
    if repeat:
        vprint("Matching repeat:\n"+repeat)

        # save output
    if saveout:
        print("\n----------------------------------------------------------------------")
        print("Saving output in: \n"+outdir)
        print("----------------------------------------------------------------------\n")
        if not os.path.exists(outdir+seqid):
                os.makedirs(outdir+seqid)
        print("Writing metadata to info_"+seqid+".txt")
        f=open(outdir+seqid+"/info_"+seqid+".txt","w")
        f.write(info_text)
        f.close()
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
        f=open(outdir+seqid+"/CRISPR_"+seqid+".ffn", "w")
        f.close()
        f=open(outdir+seqid+"/CRISPR_"+seqid+".ffn", "a")
        f.write("@made by L-F-S\n@ University of Trento\n@lorenzo.signorini@alumni.unitn.it\n\n")
        f.write("CRISPR array sequence for "+ seqid+" "+feature+".\n"+CRISPRstrandstring+"\n\n\n")
        f.write(CRISPRheader)
        f.write(CRISPRarraysequence+"\n\n")
        if repeat:
            f.write("Matching repeat:\t"+repeat)
        f.close()

        print("Saving original genome")
        full_genome_full_path, annotation_folder=originalpath.print_path(caslocus.genomename)
        old_genome_name=full_genome_full_path.split("/")[-1].rstrip(".fa")
        print(old_genome_name)
        os.popen("cp "+full_genome_full_path+" "+outdir+seqid+"/"+caslocus.genomename+".ffn") #change the original  extension .fa to .ffn because  .faa is used to dentoe aminoacidic Fasta WARNING: genomename.ffn is also the name of prokka annotation file, so those 2 should be ahndled in separat efolders (prokka annotation is temporarily processed in wdir, in fact

        print("Writing "+feature+" nucleotidic sequence to "+feature+"_seqid"+".ffn")
        # v4 21/04/2020
        for record in SeqIO.parse(annotation_folder+caslocus.orig_genomename+".ffn","fasta"):
            if record.id.startswith(seqid):
                print(record.id)
                print("++++++++++++++++++++++++++\n",record.seq)
                record.description=feature+" len="+str(len(record.seq))+" genome="+caslocus.genomename+" SGB="+str(caslocus.SGB)+" contig="+caslocus.contigname
                SeqIO.write(record, outdir+seqid+"/"+feature+"_"+seqid+".ffn","fasta")
                break

        print("Aligning Cas9 amino acid sequence with references")
        ref_fasta=wdir+"control/uniprot_working_Cas9s.fasta"
        cas9_aa_path=outdir+seqid+"/"+feature+"_"+seqid+".faa"
        alignments=[]
        for row in list(SeqIO.parse(cas9_aa_path,'fasta'))+list(SeqIO.parse(ref_fasta, 'fasta')):
            tempseq=SeqRecord(row.seq, id=row.id, description="")
            alignments.append(tempseq)

        SeqIO.write(alignments, outdir+seqid+"/msa.faa", "fasta")
        cline= ClustalwCommandline("clustalw", infile=outdir+seqid+"/msa.faa", outfile=outdir+seqid+"/msa.aln")
        os.system(str(cline))

        print("Printing SVG file of alignment from Jalview..")
        os.chdir(outdir+seqid+"/")
        os.system(" jalview -open "+outdir+seqid+"/msa.aln -nodisplay -colour Clustal -features "+wdir+\
                  "control/ref_features_colore  -svg msa.svg")






if __name__=="__main__":
    main_descr="Print information about a given effector Cas Sequence ID."
    plus="+"*len(main_descr)
    parser=argparse.ArgumentParser(description="+"*5+"\t\t"+main_descr+"\t\t"+"+"*5 )
    parser.add_argument("ID", type=str, help="sequence ID")
    parser.add_argument("-v", action="store_true", help="verbose output")
    parser.add_argument("-f", type=str, help="effector Cas name (default= Cas9)", default="Cas9")
    parser.add_argument("-o", type=str, help="output directory, default =/shares/CIBIO-Storage/CM/news/users/lorenzo.signorini/9output/<feature>/", default="/shares/CIBIO-Storage/CM/news/users/lorenzo.signorini/9output/")
    parser.add_argument("-w", type=str, help="working directory, default =/shares/CIBIO-Storage/CM/news/users/lorenzo.signorini/", default="/shares/CIBIO-Storage/CM/news/users/lorenzo.signorini/")
    parser.add_argument("-d", action="store_true",  help="Save sequence output")
    parser.add_argument("-t", type=str, help="tracRNA sequence")
    parser.add_argument("-s", type=str, help="tracrRNA strand. Possible valuesi ( +1,-1)")
    parser.add_argument("-c", type=int, help="actual CRISPRarray strand (as discovered via previous tracrRNA analysis)\
                        . Input \'+1\' if CRISPRarray strand is the same annotated from minced algorithm (forward strand), \'-1\' if reverse strand.")
    parser.add_argument("-r", type=str, help="repeat sequence matching tracrRNA")
    args=parser.parse_args()
    seqid =args.ID
    feature=args.f
    outdir=args.o+feature+"/"
    wdir=args.w
    tracrRNA=args.t
    repeat=args.r

    if not os.path.exists(outdir):
            os.makedirs(outdir)

    get_ID_info(seqid, feature,args.v,args.d,outdir,tracrRNA, repeat,wdir)
