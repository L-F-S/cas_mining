from numpy import average
from numpy import round

class CRISPR_array:
    """ A class representing a CRISPR array.
        It holds the sequences of repeats and spacers and their length.
    """

    def __init__(self, coordinates, positions, repeats, spacers, contig):
        self._start = coordinates[0]
        self._end = coordinates[1]
        self._positions = positions
        self._repeats = repeats
        self._spacers = spacers
        self._contig = contig

    def get_start(self):
        return self._start

    def get_end(self):
        return self._end

    def get_coord(self):
        return [self._start, self._end]

    def get_positions(self):
        return self._positions

    def get_repeats(self):
        return self._repeats

    def get_spacers(self):
        return self._spacers

    def get_contig(self):
        return self._contig

    def average_repeat_len(self):
        return int(round(average([len(seq) for seq in self._repeats])))

    def average_spacer_len(self):
        return int(round(average([len(seq) for seq in self._spacers])))

    def size(self):
        return self._end - self._start

def read_minced_old(filename):
    loci = []
    with open(filename, "r") as fin:
        contig = None

        coord = []
        positions = []
        repeats = []
        spacers = []

        for line in fin:
            if line.startswith("Sequence"):
                contig = line.split(' ')[1].strip('\'')
            elif line[0] == 'C':
                if len(coord) is not 0:
                    loci.append(CRISPR_array(coord, positions, repeats, spacers, contig))
                    coord = []
                    positions = []
                    repeats = []
                    spacers = []
                
                line = line.strip().split(' ')
                coord = (int(line[5]), int(line[7]))
            elif line[0].isdigit():
                data = line.split('\t')
                positions.append(data[0])
                repeats.append(data[2])
                if len(data) == 5:
                    spacers.append(data[3])
            else:
                continue
        loci.append(CRISPR_array(coord, positions, repeats, spacers, contig))
    return loci

def read_minced(filename):
    loci = []
    with open(filename, "r") as fin:
        contig = None

        coord = []
        positions = []
        repeats = []
        spacers = []

        for line in fin:
            if line.startswith("Sequence") or line[0] == 'C':
                if len(coord) != 0:
                    loci.append(CRISPR_array(coord, positions, repeats, spacers, contig))
                    coord = []
                    positions = []
                    repeats = []
                    spacers = []
                if line.startswith("Sequence"):
                    contig = line.split(' ')[1].strip('\'')
                else:
                    line = line.strip().split(' ')
                    coord = (int(line[5]), int(line[7]))
            elif line[0].isdigit():
                data = line.split('\t')
                positions.append(data[0])
                repeats.append(data[2])
                if len(data) == 5:
                    spacers.append(data[3])
            else:
                continue
        loci.append(CRISPR_array(coord, positions, repeats, spacers, contig))
    return loci

def get_repeats(repeats, rep):
    sequences = []
    for repeat in repeats:
        if repeat == "."*len(rep):
            sequences.append(rep)
        else:
            raise Exception("Unknown repeat pattern: " + repeat + " Update get_repeats!")
    return sequences

def read_pilercr(filename):
    loci = []
    with open(filename, "r") as fin:
        contig = None

        coord = []
        positions = []
        repeats = []
        spacers = []

        for line in fin:
            if line.startswith(">"):
                contig = line.strip().strip(">")
            elif len(line.strip()) > 0 and line.strip()[0].isdigit():
                data = [i for i in line.strip().split(" ") if i != ""]
                if len(data) >= 6:
                    positions.append(int(data[0]))
                    if len(data) == 7:
                        repeats.append(data[5])
                        spacers.append(data[6])
                    else:
                        repeats.append(data[4])
                elif len(data) == 4:
                    rep = data[3]
                    repeats = get_repeats(repeats, rep)
                    start = min(positions)
                    end = max(positions) + len(repeats[-1])
                    coord = (start, end)
                    loci.append(CRISPR_array(coord, positions, repeats, spacers, contig))
                    positions = []
                    repeats = []
                    spacers = []
            elif line.startswith("SUMMARY"):
                break
            else:
                continue
    return loci

class ARNold_terminator:
    """ A class to represent the output of ARNold.
    """
    def __init__(self, position, strand, method, sequence, fold, energy):
        self._position = position
        self._strand = strand
        self._method = method
        self._sequence = sequence
        self._fold = fold
        self._energy = energy

    def get_strand(self):
        return self._strand
    
    def get_method(self):
        return self._method
    
    def get_sequence(self):
        return self._sequence
    
    def get_fold(self):
        return self._fold
    
    def get_energy(self):
        return self._energy
    
    def get_position(self):
        return self._position

    def get_start(self):
        return self._position

    def get_end(self):
        return self._position + len(self._sequence)

    def __str__(self):
        return " ".join([str(self._position), self._strand, self._sequence, self._fold, str(self._energy), self._method])

def read_ARNold(filename):
    with open(filename, "r") as fin:
        fin.readline()
        terms = []
        for line in fin:
            line = [string for string in line.strip().split(' ') if string is not '']
            new = ARNold_terminator(int(line[0]), line[2], line[1], line[3].replace('(','').replace(')','').upper(),
                             line[3], line[4] if line[4] != 'NA' else None)
            terms.append(new)
        return terms

class blastn_match:
    """ A class that represents the output of BLAST used to search for antirepeats.
    """

    def __init__(self, qseqid, seqid, pident, qlen, length, mismatch, gapopen, qseq, sseq, sstart, send, evalue, sstrand):
        self._qseqid = qseqid
        self._seqid = seqid
        self._pident = pident
        self._qlen = qlen
        self._length = length
        self._mismatch = mismatch
        self._gapopen = gapopen
        self._qseq = qseq
        self._sseq = sseq
        self._sstart = sstart
        self._send = send
        self._evalue = evalue
        self._sstrand = sstrand

    def get_query_id(self):
        return self._qseqid

    def get_subject_id(self):
        return self._seqid

    def get_subject_seq(self):
        return self._sseq

    def get_subject_start(self):
        return self._sstart

    def get_subject_end(self):
        return self._send

    def get_subject_strand(self):
        return self._sstrand

    def get_subject_coord(self):
        return self._sstart, self._send

    def get_evalue(self):
        return self._evalue

    def get_subject_len(self):
        return self._send - self._sstart + 1

    # ADD MORE ACCESSORS

    def __str__(self):
        return " ".join([self._qseqid, self._seqid, str(self._pident), str(self._qlen), str(self._length), str(self._mismatch), str(self._gapopen), self._qseq, self._sseq, str(self._sstart), str(self._send), str(self._evalue), self._sstrand])

    def is_near(self, blastn_obj, r):
        coord = blastn_obj.get_subject_coord()
        c = self.get_subject_coord()
        if abs(c[0] - coord[0]) < r or abs(c[1] - coord[1]) < r:
            return True
        else:
            return False

def read_balstout(filename):
    with open(filename, "r") as fin:
        matches = []
        for line in fin:
            l = line.strip().split('\t')
            sstart = int(l[9])
            send = int(l[10])
            if sstart > send:
                sstart, send = send, sstart
            m = blastn_match(l[0], l[1], float(l[2]), int(l[3]), int(l[4]), float(l[5]), float(l[6]), l[7], l[8],
                             sstart, send, float(l[11]), "+" if l[12] == "plus" else "-")
            matches.append(m)
        return matches

def filter_blastout(blastout, dist):
    """ Filter blast ouptut removing the matches that are basically the same
        sequence.
    """
    blastout = {i:blastout[i] for i in range(len(blastout))}
    unique = []
    while blastout != {}:
        match = blastout[list(blastout.keys())[0]]
        near = [key for key in blastout if match.is_near(blastout[key], dist)]
        best = match
        for k in near:
            if blastout[k].get_evalue() < best.get_evalue():
                best = blastout[k]
            del blastout[k]
        unique.append(best)
    return unique

class RNIE_terminator:
    """ A class to represent the output of ARNold.
    """
    def __init__(self, ID, start, end, strand, sequence, bit_score):
        self._ID = ID
        self._start = start
        self._end = end
        self._strand = strand
        self._sequence = sequence
        self._bit_score = bit_score
        self._method = "RNIE"

    def get_bit_score(self):
        return self._bit_score

    def get_ID(self):
        return self._ID

    def get_strand(self):
        return self._strand
    
    def get_sequence(self):
        return self._sequence

    def get_start(self):
        return self._start

    def get_end(self):
        return self._end
    
    def get_method(self):
        return self._method

    def __len__(self):
        return len(self._sequence)

    def __str__(self):
        return " ".join([self._ID, str(self._start), str(self._end), self._strand, self._sequence, str(self._bit_score)])

def read_RNIE_gff(filename, contig):
    terminators = []
    with open(filename, "r") as fin:
        for line in fin:
            l = line.strip().split('\t')
            flank_start = int(l[0].split("_")[2])
            seq_start = flank_start + int(l[3]) - 1
            seq_end = flank_start + int(l[4])
            seq = str(contig[seq_start:seq_end].seq)
            bit_score = float(l[5])
            strand = l[6]
            new = RNIE_terminator(l[0], seq_start+1, seq_end, strand, seq, bit_score)
            terminators.append(new)
    return terminators

class tracrRNA:
    """ A class to represent a tracrRNA.
    """

    def __init__(self, seq, start, end, strand, method, repeat, crispr_strand):
        self._seq = seq
        self._start = start
        self._end = end
        self._strand = strand
        self._method = method
        self._repeat = repeat
        self._crispr_strand = crispr_strand

    def get_seq(self):
        return self._seq

    def get_start(self):
        return self._start

    def get_end(self):
        return self._end

    def get_stand(self):
        return self._strand

    def get_method(self):
        return self._method

    def get_coord(self):
        return self._start, self._end

    def get_repeat(self):
        return self._repeat

    def get_crispr_strand(self):
        return self._crispr_strand

    def __str__(self):
        return " ".join([str(self._start), str(self._end), self._strand, self._seq])

class hybrid:
    """A class to represent a tracrRNA-crRNA hybrid
    """
    def __init__(self, mfe, structure):
        self._mfe = mfe
        self._structure = structure

    def get_mfe(self):
        return self._mfe

    def get_structure(self):
        return self._structure

    def get_start(self):
        return int(self._structure[0].strip().split("  ")[1])

    def get_end(self):
        seq = []
        f = self._structure
        for i in range(11, len(f[1])):
            if f[1][i] == " ":
                if f[2][i] == " ":
                    continue
                else:
                    seq.append(f[2][i])
            elif f[1][i] == "3":
                break
            else:
                seq.append(f[1][i])
        if f[1][10] == " ":
            return len(seq)
        else:
            return len(seq) + self.get_start() - 1

def read_hybrid(filename):
    f = open(filename, "r")
    lines = f.readlines()
    f.close()
    mfe = float(lines[5].strip().split(" ")[1])
    structure = lines[8:13]
    structure[1] = "tracr "+structure[1][6:len(structure[1])]
    structure[4] = "repeat"+structure[4][6:len(structure[4])]
    return hybrid(mfe, structure)

class RNA_fold:
    """A class to represent the secondary structure of an RNA
    """

    def __init__(self, sequence, fold, mfe, method):
        self._sequence = sequence
        self._fold = fold
        self._mfe = mfe
        self._method = method

    def get_sequence(self):
        return self._sequence

    def get_fold(self):
        return self._fold

    def get_mfe(self):
        return self._mfe

    def get_method(self):
        return self._method

def read_RNAfold(filename):
    f = open(filename, "r")
    lines = f.readlines()
    f.close()
    seq = lines[1].strip()
    fold = lines[2].strip().split(" ")[0]
    mfe = float(lines[2].strip().split(" ")[-1].strip("(").strip(")"))
    return RNA_fold(seq, fold, mfe, "RNAfold")

class mfold:
    """A class to represent the output of mfold
    """
    def __init__(self, structure, mfe):
        self._structure = structure
        self._mfe = mfe
    
    def get_structure(self):
        return self._structure
    
    def get_mfe(self):
        return self._mfe
    
    def __str__(self):
        return "mfe = {} kcal/mol\n\n".format(self._mfe)+"".join(self._structure)

def read_mfold(filename):
    structures = []
    with open(filename, "r") as fin:
        structure = []
        mfe = None
        for line in fin:
            if len(line.strip()) != 0:
                if line.strip().startswith("Structure"):
                    if int(line.strip().split(" ")[-1]) > 1:
                        structures.append(mfold(structure, mfe))
                        structure = []
                elif line.strip().startswith("Initial"):
                    mfe = float(line.strip().split(" ")[-1])
                elif line.strip().startswith("Folding"):
                    continue
                else:
                    structure.append(line)
        structures.append(mfold(structure, mfe))
    return structures
