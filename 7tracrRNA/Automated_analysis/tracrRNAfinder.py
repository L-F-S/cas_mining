#!/usr/bin/env python3

import os
import sys
import argparse
import re
from tempfile import TemporaryDirectory
from numpy import unique
from pandas import read_csv
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Alphabet import generic_dna, generic_protein, RNAAlphabet
from Bio import pairwise2
from CRISPR_class import *

FEATURE="Cas9"
LOCUS_DIST=3500
TERMINATOR_DIST=200
MIN_DIST_FROM_CRISPR=50
MIN_POLYT_DIST=40
MIN_POLYT_LEN=3

def params():
    parser = argparse.ArgumentParser(description="Searches for one or more tracrRNAs near a CRISPR locus, specified"
                                                 " by the ID of an effector protein (a Cas9 by default).",
                                     epilog="Author: Matteo Ciciani")
    parser.add_argument("effector_ID", type=str, help="The ID of the effector protein")
    parser.add_argument("genome_directory", type=str, help="The path to the directory containing the genome")
    parser.add_argument("effector_dataset", type=str, help="The path to the effector variants table")
    parser.add_argument("-f", "--feature", type=str, required=False,
                        help="The type of effector protein. Default is Cas9", default=FEATURE)
    parser.add_argument("-o", "--output-directory", type=str, required=False,
                        help="The path to the output directory. Defaults to the current directory")
    parser.add_argument("-i", "--intermediate-results", action="store_true",
                        help="Save intermediate results in the output directory")
    parser.add_argument("-m", "--minced-directory", type=str, required=False,
                        help="The path to the directory containing the minced output file. Defaults to the genome"
                             " directory")
    parser.add_argument("-l", "--locus-dist", type=int, required=False,
                        help="The max distance between the locus boundaries and the anti-repeat"
                             " Defaults to {}".format(LOCUS_DIST), default=LOCUS_DIST)
    parser.add_argument("-t", "--terminator-dist", type=int, required=False,
                        help="The max distance between the anti-repeat and the terminator"
                             " Defaults to {}".format(TERMINATOR_DIST), default=TERMINATOR_DIST)
    parser.add_argument("-v", "--verbose", action="store_true", help="Be verbose")
    parser.add_argument("-g", "--no-graph", action="store_true", help="Suppress graphical output")
    parser.add_argument("-s", "--on-server", action="store_true", help="Set if run on server")
    parser.add_argument("--min-crispr-dist", type=int, required=False,
                        help="The min distance between the CRISPR array and the anti-repeat"
                             " Defaults to {}".format(MIN_DIST_FROM_CRISPR), default=MIN_DIST_FROM_CRISPR)
    parser.add_argument("--min-polyT-len", type=int, required=False,
                        help="The min length of the poly-T at the end of a terminator."
                             " Defaults to {}".format(MIN_POLYT_LEN), default=MIN_POLYT_LEN)
    parser.add_argument("--min-polyT-dist", type=int, required=False,
                        help="The min distance between an anti-repeat and a poly-T."
                             " Defaults to {}".format(MIN_POLYT_DIST), default=MIN_POLYT_DIST)
    parser.add_argument("--allow-internal", action="store_true", help="Allow blastn matches inside the "
                                                                      "CRISPR array")
    parser.add_argument("--all-analyses", action="store_true", help="Search tracrRNAs using every method")
    parser.add_argument("--improve-loop", action="store_true", help="Attempts to improve tracrRNAs with a "
                                                "bad loop in the upper stem of the repeat:anti-repeat stem loop")
    args = parser.parse_args()
    return args

def check_params(pars):
    if os.path.exists(pars.genome_directory):
        if not os.path.isdir(pars.genome_directory):
            sys.stderr.write("Error: {} is not a directory!\n".format(pars.genome_directory))
            sys.exit(1)
    else:
        sys.stderr.write("Error: {} does not exist!\n".format(pars.genome_directory))
        sys.exit(1)
    if os.path.exists(pars.effector_dataset):
        if not os.path.isfile(pars.effector_dataset):
            sys.stderr.write("Error: {} is not a file!\n".format(pars.effector_dataset))
            sys.exit(1)
    else:
        sys.stderr.write("Error: {} does not exist!\n".format(pars.effector_dataset))
        sys.exit(1)
    if pars.output_directory is not None:
        if os.path.exists(pars.output_directory):
            if not os.path.isdir(pars.output_directory):
                sys.stderr.write("Error: {} is not a directory!\n".format(pars.output_directory))
                sys.exit(1)
    if pars.minced_directory is not None:
        if os.path.exists(pars.minced_directory):
            if not os.path.isdir(pars.minced_directory):
                sys.stderr.write("Error: {} is not a directory!\n".format(pars.minced_directory))
                sys.exit(1)
        else:
            sys.stderr.write("Error: {} does not exist!\n".format(pars.minced_directory))
            sys.exit(1)
    if pars.locus_dist <= 0:
        sys.stderr.write("Warning: locus distance must be positive, using default value {}\n".format(LOCUS_DIST))
        pars.locus_dist = LOCUS_DIST
    if pars.terminator_dist <= 0:
        sys.stderr.write("Warning: terminator distance must be positive, using default "
                         "value {}\n".format(TERMINATOR_DIST))
        pars.terminator_dist = TERMINATOR_DIST
    if pars.min_crispr_dist <= 0:
        sys.stderr.write("Warning: min distance from CRISPR array must be positive, using default "
                         "value {}\n".format(MIN_DIST_FROM_CRISPR))
        pars.min_crispr_dist = MIN_DIST_FROM_CRISPR
    if pars.min_polyT_len <= 0:
        sys.stderr.write("Warning: min poly-T length must be positive, using default "
                         "value {}\n".format(MIN_POLYT_LEN))
        pars.min_polyT_len = MIN_POLYT_LEN
    if pars.min_polyT_dist <= 0:
        sys.stderr.write("Warning: min poly-T distance must be positive, using default "
                         "value {}\n".format(MIN_POLYT_DIST))
        pars.min_polyT_dist = MIN_POLYT_DIST

def get_actual_outdir(outdir, feature, seqid):
    if outdir is None:
        outdir = "./"
    if outdir[-1] != "/":
        outdir = outdir+"/"
    outdir = outdir+feature+"/"+seqid+"/"
    if not os.path.exists(outdir):
        os.makedirs(outdir)
    return outdir

def get_genome_filename(genomedir, genomename, onserver, sgb):
    if genomedir[-1] != "/":
        genomedir = genomedir+"/"
    if onserver:
        genomedir = genomedir+str(sgb)+"/"
    genomefilename = genomedir+genomename+".fa"
    if os.path.exists(genomefilename):
        if not os.path.isfile(genomefilename):
            sys.stderr.write("Error: {}".format(genomefilename)+" is not a file!\n")
            sys.exit(1)
    else:
        sys.stderr.write("Error: {}".format(genomefilename) + " does not exixt!\n")
        sys.exit(1)
    return genomefilename

def get_minced_filename(mincedir, genomename, onserver, dataset):
    if mincedir[-1] != "/":
        mincedir = mincedir+"/"
    if onserver:
        mincedir = mincedir+dataset+"/"
    mincedfilename = mincedir+genomename+".fa.minced.out"
    if os.path.exists(mincedfilename):
        if not os.path.isfile(mincedfilename):
            sys.stderr.write("Error: {}".format(mincedfilename)+" is not a file!\n")
            sys.exit(1)
    else:
        sys.stderr.write("Error: {}".format(mincedfilename) + " does not exixt!\n")
        sys.exit(1)
    return mincedfilename

def run_blastn(unique_repeats, contig, genomename, inter, outdir):
    if inter:
        blast_folder = outdir + "intermediate/blastn/"
        if not os.path.exists(blast_folder):
            os.makedirs(blast_folder)
    else:
        blast_folder_obj = TemporaryDirectory()
        blast_folder = blast_folder_obj.name + "/"

    # save unique repeats in fasta file
    repeatfile = open(blast_folder+"repeat_seq.fasta", "w")
    for n, repeat in enumerate([Seq(sequence) for sequence in unique_repeats]):
        repeatfile.write(">cr_rpt_" + str(n) + "|" + genomename + "\n" + str(repeat) + "\n")
    repeatfile.close()

    # write contig to file
    contigfilename = blast_folder + "contig.fasta"
    SeqIO.write(contig, contigfilename, "fasta")

    # create blast db
    os.system("makeblastdb -in " + contigfilename + " -parse_seqids  -dbtype nucl")

    # perform blastn
    blastoutfile = blast_folder + genomename + ".trcr.blastout"
    query = blast_folder + "repeat_seq.fasta"
    blastn_command = "blastn -out " + blastoutfile + " -outfmt \"6  qseqid sseqid pident qlen length mismatch" \
                                                     " gapopen qseq sseq sstart send evalue sstrand\" -query " + \
                     query + " -db " + contigfilename + " -evalue 0.1 -word_size 11 -task blastn-short -gapopen 2" \
                                                        " -gapextend 1 -penalty -1 -reward 1"
    os.system(blastn_command)

    # read blastout
    matches = read_balstout(blastoutfile)

    # cleanup
    if not inter:
        blast_folder_obj.cleanup()

    return matches

def is_near_locus(match, cas_position, cutoff, min_dist, internal):
    locus_start = min([pos[0] for pos in cas_position.values()])
    locus_end = max([pos[1] for pos in cas_position.values()])
    CRISPR_range = range(cas_position['CRISPR'][0]-min_dist, cas_position['CRISPR'][1]+min_dist+1)
    if locus_start - cutoff < match.get_subject_start() and match.get_subject_end() < locus_end + cutoff:
        if (match.get_subject_start() not in CRISPR_range and match.get_subject_end() not in CRISPR_range) or internal:
            return True
    return False

def run_RNIE(flanks, inter, outdir, contig):
    if inter:
        RNIE_folder = outdir + "intermediate/"
        if not os.path.exists(RNIE_folder):
            os.makedirs(RNIE_folder)
    else:
        RNIE_folder_obj = TemporaryDirectory()
        RNIE_folder = RNIE_folder_obj.name + "/"
    current_dir = os.getcwd()
    os.chdir(RNIE_folder)

    #write flanks to file
    with open("flanks.fasta", "w") as fout:
        for f in flanks:
            fout.write(f + '\n' + str(flanks[f].seq) + '\n')

    # run RNIE
    os.system("$RNIE/rnie.pl --gene -f flanks.fasta")

    # extract sequences form contig
    terminators = read_RNIE_gff("flanks-geneMode-rnie.gff", contig)

    #cleanup
    os.chdir(current_dir)
    if not inter:
        RNIE_folder_obj.cleanup()

    return terminators

def find_pairs(selected, terminators):
    pairs = []
    for ID in selected:
        for t in terminators:
            if t.get_ID().split("_")[0] == ID:
                if (t.get_ID().split("_")[1] == "up" and t.get_strand() == "-") or \
                        (t.get_ID().split("_")[1] == "down" and t.get_strand() == "+"):
                    pairs.append([selected[ID], t])
    return pairs

def get_tracrRNAs(pairs, contig, unique_repeats, improve_loop):
    tracrRNAs = []
    for c in pairs:
        antirep = c[0]
        term = c[1]

        # infer crispr strand
        if antirep.get_subject_strand() == term.get_strand():
            crispr_strand = '-'
        else:
            crispr_strand = '+'

        # get actual repeat
        rep = unique_repeats[int(antirep.get_query_id().split("|")[0].split("_")[2])]
        if crispr_strand == '-':
            rep = str(Seq(rep).reverse_complement())

        # attempts to improve upper stem loop
        offset = 0
        if improve_loop and len(rep) - antirep.get_subject_len() > 3:
            r = Seq(rep).reverse_complement()
            ar = contig[antirep.get_subject_start() - 1:antirep.get_subject_end()].seq
            if term.get_strand() == "-":
                ar = ar.reverse_complement()
            al = pairwise2.align.localms(ar, r, 1, -1, -2, -1)
            al = pairwise2.format_alignment(*al[0])
            i = 0
            while al[i] == "-": i += 1
            if i > 3:
                offset = i

        # get tracrRNA sequence
        if term.get_strand() == '+':
            start = antirep.get_subject_start() - offset
            end = term.get_end()
            seq = str(contig[start - 1:end].seq)
        else:
            start = term.get_start()
            end = antirep.get_subject_end() + offset
            seq = str(contig[start - 1:end].reverse_complement().seq)

        new = tracrRNA(seq, start, end, term.get_strand(), term.get_method(), rep, crispr_strand)
        tracrRNAs.append(new)
    return tracrRNAs

def count_stem_loops(struct):
    n = 0
    previous = 0
    opened = 0
    for char in struct:
        if char == "(":
            opened += 1
        if char == ")":
            opened -= 1
        if opened == 0 and previous == 1:
            n += 1
        previous = opened
    return n

def find_polyT(seq, min_polyT_dist, min_polyT_len):
    pattern = re.compile("T"*min_polyT_len + "[T]*")
    ends = []
    i = min_polyT_dist
    m = pattern.search(seq, i)
    while m is not None:
        ends.append(m.end())
        m = pattern.search(seq, m.end())
    return ends

def analyze_fold(tracrRNAs, inter, outdir, nograph):
    if inter:
        fold_folder = outdir + "intermediate/"
        if not os.path.exists(fold_folder):
            os.makedirs(fold_folder)
    else:
        fold_folder_obj = TemporaryDirectory()
        fold_folder = fold_folder_obj.name + "/"
    mfold_folder = fold_folder + "mfold/"
    if not os.path.exists(mfold_folder):
        os.makedirs(mfold_folder)
    current_dir = os.getcwd()
    os.chdir(fold_folder)

    folds = []

    for i in range(len(tracrRNAs)):
        # get repeat and tracrRNA
        tracr = tracrRNAs[i].get_seq()
        rep = tracrRNAs[i].get_repeat()

        # use RNAhybrid to predict repeat-antirepeat base pairing
        if nograph:
            hybrid_command = "RNAhybrid " + tracr + " " + rep + " -s 3utr_human > hybrid_" + str(i + 1) + ".txt"
        else:
            hybrid_command = "RNAhybrid " + tracr + " " + rep + " -s 3utr_human -g ps > hybrid_" + str(i + 1) + ".txt"
        os.system(hybrid_command)
        hyb = read_hybrid("hybrid_" + str(i + 1) + ".txt")

        # use RNAfold to fold the tail
        tail = tracr[hyb.get_end() - 1:len(tracr)]
        tailfile = "input_tail_" + str(i + 1) + ".fasta"
        with open(tailfile, "w") as fout:
            fout.write(">tail_" + str(i + 1) + "\n" + tail)
        if nograph:
            fold_command = "RNAfold " + tailfile + " --noLP --noPS > output_tail_" + str(i + 1) + ".txt"
        else:
            fold_command = "RNAfold " + tailfile + " --noLP > output_tail_" + str(i + 1) + ".txt"
        os.system(fold_command)
        tail_fold = read_RNAfold("output_tail_" + str(i + 1) + ".txt")
        if not nograph:
            os.system("mv command_line_command_line_1.ps hybrid_" + str(i + 1) + ".ps")

        # run mfold to analyze sgRNA
        os.chdir("mfold")
        sgRNA = str(Seq(rep + "GAAA" + tracr).transcribe())
        sgRNA_file = "sgRNA_" + str(i + 1) + ".fasta"
        with open(sgRNA_file, "w") as fout:
            fout.write(">sgRNA_" + str(i + 1) + "\n" + sgRNA)
        os.system("mfold SEQ=" + sgRNA_file)
        sgRNA_folds = read_mfold("sgRNA_" + str(i + 1) + ".out")
        os.system("mv sgRNA_" + str(i + 1) + "_1.ps " + "..")
        os.chdir("..")

        folds.append([hyb, tail_fold, sgRNA_folds[0]])

    # cleanup
    os.chdir(current_dir)
    if not nograph:
        os.system("mv "+fold_folder+"*.ps " + outdir)

    if not inter:
        fold_folder_obj.cleanup()

    return folds

def has_correct_structure(vienna):
    n = 0
    opened = 0
    previous = 0
    closing = False
    for char in vienna:
        if char == "(":
            opened += 1
            if closing:
                return False
        if char == ")":
            opened -= 1
            closing = True
        if opened == 0 and previous == 1:
            n += 1
            closing = False
        previous = opened
    if n >= 2:
        return True
    return False

def has_valid_fold(tracrRNA, inter, outdir, n):
    if inter:
        fold_folder = outdir + "intermediate/flank_analysis/sgRNA_{}/".format(n)
        if not os.path.exists(fold_folder):
            os.makedirs(fold_folder)
    else:
        fold_folder_obj = TemporaryDirectory()
        fold_folder = fold_folder_obj.name + "/"
    if not os.path.exists(fold_folder):
        os.makedirs(fold_folder)
    current_dir = os.getcwd()
    os.chdir(fold_folder)

    tracr = tracrRNA.get_seq()
    rep = tracrRNA.get_repeat()
    sgRNA = str(Seq(rep + "GAAA" + tracr).transcribe())
    sgRNA_file = "sgRNA_{}.fasta".format(n)
    with open(sgRNA_file, "w") as fout:
        fout.write(">sgRNA_{}\n".format(n) + sgRNA)
    os.system("mfold SEQ=" + sgRNA_file)
    os.system("ct2dot sgRNA_{}_1.ct -1 sgRNA_{}_1.txt".format(n, n))
    fold = [line.strip() for line in open("sgRNA_{}_1.txt".format(n), "r")][2]

    # cleanup
    os.chdir(current_dir)
    if not inter:
        fold_folder_obj.cleanup()

    return has_correct_structure(fold)

def analyze_flanking_regions(selected, flanks, min_polyT_dist, min_polyT_len, contig, unique_repeats, inter, outdir, improve_loop):
    selected_tracrs = []
    n = 0
    for f in flanks:
        # get selected match
        antirep = selected[f.split("_")[0].strip(">")]

        # get correctly oriented flank sequence
        strand = "-" if f.split("_")[1] == "up" else "+"
        if strand == "-":
            flanks[f] = flanks[f].reverse_complement()

        # look for polyT in flank
        ends = find_polyT(str(flanks[f].seq), min_polyT_dist, min_polyT_len)

        # infer CRISPR strand and get correct repeat
        CRISPR_strand = '-' if antirep.get_subject_strand() == strand else '+'
        rep = unique_repeats[int(antirep.get_query_id().split("|")[0].split("_")[2])]
        if CRISPR_strand == '-':
            rep = str(Seq(rep).reverse_complement())

        # attempts to improve upper stem loop
        offset = 0
        if improve_loop and len(rep) - antirep.get_subject_len() > 3:
            r = Seq(rep).reverse_complement()
            ar = contig[antirep.get_subject_start() - 1:antirep.get_subject_end()].seq
            if strand == "-":
                ar = ar.reverse_complement()
            al = pairwise2.align.localms(ar, r, 1, -1, -2, -1)
            al = pairwise2.format_alignment(*al[0])
            i = 0
            while al[i] == "-": i += 1
            if i > 3:
                offset = i

        # extract tracr sequences
        putative_tracrs = []
        for end in ends:
            # get tracrRNA sequence
            if strand == "+":
                t_start = antirep.get_subject_start() - 1 - offset
                t_end = antirep.get_subject_end() + end
                t_seq = str(contig[t_start:t_end].seq)
            else:
                t_start = antirep.get_subject_start() - 1 - end
                t_end = antirep.get_subject_end() + offset
                t_seq = str(contig[t_start:t_end].reverse_complement().seq)
            putative_tracrs.append(tracrRNA(t_seq, t_start, t_end, strand, "Secondary structure analysis", rep,
                                            CRISPR_strand))

        # predict secondary structure of each tracrRNA
        for tracr in putative_tracrs:
            if has_valid_fold(tracr, inter, outdir, n):
                selected_tracrs.append(tracr)
            n += 1

    return selected_tracrs

def write_output(outdir, seqid, feature, cas_dataset, cas_position, f_strand, prot_seq, tracrRNAs, folds):
    fout = open(outdir + "tracrRNA_" + seqid + ".txt", "w")
    fout.write("-" * 80 + "\n")
    fout.write("tracrRNA sequences for " + feature + " protein ID " + seqid + "\n")
    fout.write("-" * 80 + "\n")
    fout.write("Contig name: " + cas_dataset[cas_dataset["Seq ID"] == seqid]["Contig"].iloc[0] + "\n")
    fout.write("Dataset: " + cas_dataset[cas_dataset["Seq ID"] == seqid]["Study"].iloc[0] + "\n")
    fout.write("Bin: " + cas_dataset[cas_dataset["Seq ID"] == seqid]["Genome Name"].iloc[0] + "\n")
    fout.write("SGB: " + str(cas_dataset[cas_dataset["Seq ID"] == seqid]["SGB ID"].iloc[0]) + "\n")
    fout.write("CRISPR position: " + str(cas_position["CRISPR"][0]) + " " + str(cas_position["CRISPR"][1]) + "\n")
    fout.write("Cas2 position: " + str(cas_position["Cas2"][0]) + " " + str(cas_position["Cas2"][1]) + "\n")
    fout.write("Cas1 position: " + str(cas_position["Cas1"][0]) + " " + str(cas_position["Cas1"][1]) + "\n")
    fout.write(feature + " position: " + str(cas_position[feature][0]) + " " + str(cas_position[feature][1]) + "\n")
    fout.write(feature + " strand: " + f_strand + "\n")
    fout.write(feature + " length: " + str(len(prot_seq) - 1) + " aa" + "\n")
    fout.write("Number of tracrRNAs found: " + str(len(tracrRNAs)) + "\n")
    fout.write("-" * 80 + "\n")

    # print tracrRNAs data
    for i in range(len(tracrRNAs)):
        fout.write("-" * 80 + "\n")
        fout.write("tracrRNA_" + str(i + 1) + "\n\n")
        fout.write("Position: " + str(tracrRNAs[i].get_start()) + " " + str(tracrRNAs[i].get_end()) + "\n")
        fout.write("Strand: " + tracrRNAs[i].get_stand() + "\n")
        fout.write("Length: " + str(len(tracrRNAs[i].get_seq())) + "\n")
        fout.write("Sequence:\n")
        j = 0
        s = tracrRNAs[i].get_seq()
        while j < len(s):
            fout.write(s[j:j + 80] + "\n")
            j += 80
        fout.write("Inferred CRISPR array strand: " + tracrRNAs[i].get_crispr_strand() + "\n")
        rep = tracrRNAs[i].get_repeat()
        fout.write("Matched repeat: " + rep + "\n")
        fout.write("Method: " + tracrRNAs[i].get_method() + "\n")
        fout.write("sgRNA sequence:\n")
        sgRNA = str(Seq(rep + "GAAA" + tracrRNAs[i].get_seq()).transcribe())
        j = 0
        while j < len(sgRNA):
            fout.write(sgRNA[j:j + 80] + "\n")
            j += 80
        fout.write("\nRepeat - antirepeat hybrid structure:\n")
        for l in folds[i][0].get_structure(): fout.write(l)
        fout.write("Repeat - antirepeat hybrid mfe: " + str(folds[i][0].get_mfe()) + " kcal/mol\n")
        fout.write("\nTail secondary structure:\n")
        fout.write(folds[i][1].get_sequence() + "\n")
        fout.write(folds[i][1].get_fold() + "\n")
        fout.write("Number of stem loops: " + str(count_stem_loops(folds[i][1].get_fold())) + "\n")
        fout.write("Tail secondary structure mfe: " + str(folds[i][1].get_mfe()) + " kcal/mol\n")
        fout.write("\nsgRNA secondary structure:\n " + "".join(folds[i][2].get_structure()))
        fout.write("sgRNA secondary structure mfe: " + str(folds[i][2].get_mfe()) + " kcal/mol\n")

    # parameters
    fout.close()

def perform_analysis(pars):
    seqid = pars.effector_ID
    genomedir = pars.genome_directory
    effector_data = pars.effector_dataset
    feature = pars.feature
    outdir = pars.output_directory
    inter = pars.intermediate_results
    mincedir = genomedir if pars.minced_directory is None else pars.minced_directory
    locus_dist = pars.locus_dist
    terminator_dist = pars.terminator_dist
    verb = pars.verbose
    nograph = pars.no_graph
    onserver = pars.on_server
    min_dist = pars.min_crispr_dist
    min_polyT_dist = pars.min_polyT_dist
    min_polyT_len = pars.min_polyT_len
    internal = pars.allow_internal
    all_analyses = pars.all_analyses
    improve_loop = pars.improve_loop

    # set up output directory
    if verb:
        print("Starting analysis for " + feature + " protein ID: " + seqid + "!")
        print("Setting up output directory...")
    outdir = get_actual_outdir(outdir, feature, seqid)

    # read and check input data
    if verb:
        print("Reading input files...")
        print("Reading effector dataset...")
    cas_dataset = read_csv(effector_data, index_col=0)
    try:
        genomename = cas_dataset[cas_dataset["Seq ID"] == seqid]["Genome Name"].iloc[0]
        dataset = cas_dataset[cas_dataset["Seq ID"] == seqid]["Study"].iloc[0]
        contigname = cas_dataset[cas_dataset["Seq ID"] == seqid]["Contig"].iloc[0]
        sgb = cas_dataset[cas_dataset["Seq ID"]==seqid]["SGB ID"].iloc[0]
        cas_position = {feature: [int(n) for n in cas_dataset[cas_dataset["Seq ID"] == seqid]["Pos"].iloc[0].split()],
                        "Cas1": [int(n) for n in
                                 eval(cas_dataset[cas_dataset["Seq ID"] == seqid]["prokka_cas1"].iloc[0])[0][1:]],
                        "Cas2": [int(n) for n in
                                 eval(cas_dataset[cas_dataset["Seq ID"] == seqid]["prokka_cas2"].iloc[0])[0][1:]],
                        "CRISPR": [int(n) for n in
                                   eval(cas_dataset[cas_dataset["Seq ID"] == seqid]["minced_CRISPR"].iloc[0])[0][1:]]
                        }
    except IndexError:
        sys.stderr.write("Index error: check whether sequence ID "+seqid+" is present in {}!\n".format(effector_data))
        sys.exit(1)
    except KeyError as err:
        sys.stderr.write("Error: key {}".format(err)+" is not present in {}!\n".format(effector_data))
        sys.exit(1)

    if verb:
        print("Retrieving file names...")
    genomefilename = get_genome_filename(genomedir, genomename, onserver, sgb)
    mincedCRISPRfilename = get_minced_filename(mincedir, genomename, onserver, dataset)

    if verb:
        print("Reading minced file...")
    CRISPR_arrays = read_minced(mincedCRISPRfilename)
    if any([array.get_contig() == contigname for array in CRISPR_arrays]):
        if any([array.get_coord() == cas_position['CRISPR'] for array in CRISPR_arrays]):
            array = [l for l in CRISPR_arrays if l.get_coord() == cas_position['CRISPR'] and
                     l.get_contig() == contigname][0]
        else:
            sys.stderr.write("Error: no CRISPR array with coordinates: {}, {}!\n".format(cas_position['CRISPR'][0],
                                                                                      cas_position['CRISPR'][1]))
            sys.exit(1)
    else:
        sys.stderr.write("Error: Cas9 and CRISPR array are on different contigs!\n")
        sys.exit(1)
    unique_repeats = unique(array.get_repeats())

    if verb:
        print("Reading genome file...")
    contig = SeqIO.index(genomefilename, "fasta")[array.get_contig()]

    # get feature strand
    if verb:
        print("Checking {} strand...".format(feature))
    f_start = cas_position[feature][0]
    f_end = cas_position[feature][1]
    f_seq = contig[f_start - 1:f_end]
    prot_seq = cas_dataset[cas_dataset["Seq ID"] == seqid]["Seq"].iloc[0] + "*"
    if str(f_seq.translate().seq)[1:] == prot_seq[1:]:
        f_strand = "+"
        if str(f_seq.translate().seq)[0] != prot_seq[0]:
            sys.stderr.write("Warning: {} has unusual START codon {}!\n".format(feature, str(f_seq[0:3].seq)))
            f_strand = f_strand + " Warning: {} has unusual START codon {}!\n".format(feature, str(f_seq[0:3].seq))
    elif str(f_seq.reverse_complement().translate().seq)[1:] == prot_seq[1:]:
        f_strand = "-"
        if str(f_seq.reverse_complement().translate().seq)[0] != prot_seq[0]:
            sys.stderr.write("Warning: {} has unusual START codon {}!\n".format(feature, str(f_seq.reverse_complement()[
                                                                                       0:3].seq)))
            f_strand = f_strand + " Warning: {} has unusual START codon {}!\n".format(feature, str(
                f_seq.reverse_complement()[0:3].seq))
    else:
        sys.stderr.write("Warning: {} sequence in data table is not identical to the one in the genome!\n".format(
            feature))
        f_strand = "NA"

    # run blastn, read results, filter them and select repeats near the CRISPR locus
    if verb:
        print("Running blastn...")

    matches = run_blastn(unique_repeats, contig, genomename, inter, outdir)
    matches = [match for match in matches if match.get_subject_id() == array.get_contig()]
    if matches == []:
        sys.stderr.write("Warning: blastn found no matches in the contig!\n")
        sys.exit(0)
    filtered_matches = filter_blastout(matches, 3)
    selected = {str(i): match for i, match in enumerate(filtered_matches) if
                is_near_locus(match, cas_position, locus_dist, min_dist, internal)}
    if selected == {}:
        sys.stderr.write("Warning: no antirepeats present near the locus! Current locus-dist = {}\n".format(locus_dist))
        sys.exit(0)
    if verb:
        print("{} antirepeats selected.".format(len(selected)))

    # define flanking regions
    if verb:
        print("Defining flanking regions...")
    flanks = {}
    for ID in selected:
        start = max(selected[ID].get_subject_start() - 1 - terminator_dist, 0)
        end = selected[ID].get_subject_start() - 1
        flanks[">" + ID + "_up_" + str(start) + "_" + str(end)] = contig[start:end]
    for ID in selected:
        start = selected[ID].get_subject_end()
        end = min(selected[ID].get_subject_end() + terminator_dist, len(contig))
        flanks[">" + ID + "_down_" + str(start) + "_" + str(end)] = contig[start:end]

    # run RNIE
    if verb:
        print("Running RNIE...")
    terminators = run_RNIE(flanks, inter, outdir, contig)
    if terminators == []:
        print("Warning: no terminator found near the blastn matches! Current terminator-dist = {}\n".format(
            terminator_dist))

    # find valid antirepeat - terminator pairs
    if verb:
        print("Finding valid antirepeat - terminator pairs...")
    pairs = find_pairs(selected, terminators)
    if verb:
        if pairs == []:
            print("No valid antirepeat - terminator pair found!")

    # Extract tracrRNAs
    if verb:
        print("Extracting tracrRNAs...")
    tracrRNAs = []
    if pairs != []:
        tracrRNAs.extend(get_tracrRNAs(pairs, contig, unique_repeats, improve_loop))
    if pairs == [] or all_analyses:
        if verb:
            print("Analyzing flanking regions...")
        tracrRNAs.extend(analyze_flanking_regions(selected, flanks, min_polyT_dist, min_polyT_len, contig,
                                                unique_repeats, inter, outdir, improve_loop))
    if tracrRNAs == []:
        sys.stderr.write("Warning: no tracrRNA found near the locus!\n")
        sys.exit(0)

    # analyze tracrRNAs secondary structure
    if verb:
        print("Analyzing tracrRNAs secondary structure...")
    folds = analyze_fold(tracrRNAs, inter, outdir, nograph)

    if verb:
        print("Writing output to file...")
    write_output(outdir, seqid, feature, cas_dataset, cas_position, f_strand, prot_seq, tracrRNAs, folds)

    if verb:
        print("Done!")

if __name__ == "__main__":
    args = params()
    check_params(args)
    perform_analysis(args)
