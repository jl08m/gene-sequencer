# FILE: genseq.py
# AUTHOR: Jia Li
# CREATE DATE: 2012-07-10
# VERSION: 1.2
# CHANGE LOG: Added json dump to file

#########################################
#This can be used as a stand alone script
#if all that is needed is the full list of sequences

# Generates a list of permutated sequences based upon the parameters
# Parameters:
# bases: The number of bases to change per permutation.
# spacer length: The number of bases that remain unchanged within the bases that are changed.
# spacer position: The position within the changed bases to start the unchanged spacer.
# start position: position to start permutation in the whole sequence.
# end position: position to end permutation in the whole sequence.
#########################################
from sequence import *
import sys
from re import split as re

HEADER_FORM = 2
HEADER_SPLIT = ' |=|:|-'
BASES = ['A', 'C', 'G', 'T']
NUM_BASES = len(BASES)
_sequence = ''
_seq_length = 0

# example usage: geneseq.py somesequence.fasta 4,2,2,0,0 output.txt
USAGE = "USAGE: gen-seq.py <sequence> <bases, spacer length, spacer position, beginning, end]> \
         more help:  gen-seq.py help"
HELP = """  Sequence file format:
            sequence ID string
            sequence (may be in multiple lines with line breaks)

            Example:
            >hg19_dna range=chr11:108093980-108094130 5'pad=0 3'pad=0 strand=+ repeatMasking=none
            TCCTCTCCCCAGACCGCCAATCTCATGCACCCCTCCAGAGTGGCCCTTGA
            CTCCTCCCTCTCCTCACTCCATCTTTCCTGGCCTCTCTCCGGGTGCTTAG
            CGGACTTGGCCAATAACCTCCTCCTTTTAAACGCCCTGAATTGAACCCTG
            C

            Parameter sequence:
            number of bases, spacer length, spacer position, start position, end position
            Bases is required, put 0 for default behavior.

            number of bases --  The number of bases at a time that go through permutation. Referred to as segment.
            spacer length   --  The number of bases that a 'spacer' will have. Spacer is a section that will remain
                                unchanged in segment during permutation.
                                Default behavior is 0 length.
            spacer position --  Which position in the segment the spacer will start. Spacer length must be specified.
                                Default behavior is 1.
            start position  --  Which position in the original sequence to start the segment
                                Default position is 1.
            end_position    --  Which position in the original sequence to stop the segment.
                                Default position is the end of the sequence.

            Example:
                Parameters: 4,2,2,3,8

                Position:   1 2 3 4 5 6 7 8 9 10
                Base:       A C G T A C G T A C

                Positions that get ignored      1 2 9 10
                                                A C A C

                Loop through segments.
                Spacer position starts at 2 of segment, length of 2. Sequence positions 4,5 remain unchanged.
                                                Original        Example permutation
                Positions for 1st segment       3 4 5 6         3 4 5 6     3 4 5 6
                                                G T A C         T T A C     T T A G

                Positions for 2nd segment       4 5 6 7         4 5 6 7     4 5 6 7
                                                T A C G         A A C G     A A C C
        """


def validate_parameters(raw_parameters):
    parameters = []
    try:
        if len(raw_parameters) != 5:
            raise IndexError('Incorrect number of parameters. Expected 5, actual', len(raw_parameters))
        for raw_arg in raw_parameters:
            arg = int(raw_arg, 10)
            if arg < 0:
                raise ValueError('Argument ', raw_arg, ' cannot be negative')
            parameters.append(arg)

        bases = parameters[0]
        spacer_length = parameters[1]
        spacer_position = parameters[2]
        start = parameters[3]
        end = parameters[4]

        if bases < 1:
            raise ValueError("Bases to change cannot be less than 1")
        if bases > _seq_length:
            raise ValueError("Bases to change cannot be greater than length of sequence")

        if spacer_length > (bases-2):
            raise ValueError("Spacer length exceeds allowed bases-2")
        if spacer_position + spacer_length >= bases:
            raise ValueError("Spacer position is invalid with given spacer lenght")

        if start < 1:
            raise ValueError("Start Position cannot be less than 1")
        if start > _seq_length:
            raise ValueError("Start Position cannot be greater than length of Sequence")

        if end < start:
            raise ValueError("End Position cannot be less than start position")
        if end > _seq_length:
            raise ValueError("End Position cannot be more than length of Sequence")

    except TypeError, e:
        sys.stderr('TypeError:', e)
    except ValueError, e:
        sys.stderr('ValueError:', e)
    except IndexError, e:
        sys.stderr('IndexError:', e)

    return parameters


def get_sequence(sequence_filename):
    try:
        sequence_file = open(sequence_filename, 'r')
        sequence_id = sequence_file.readline()
        global _sequence
        _sequence = sequence_file.read().replace('\n', '')
        global _seq_length
        _seq_length = len(_sequence)
        sequence_file.close()

    except IOError, e:
        sys.stderr('Error encountered opening file:', e)

def print_header(sequence_id):
    header = sequence_id.split()
    gene_id = '>' + header[0]
    rangefield = re.split(HEADER_SPLIT, header[1])
    chrom_id = rangefield[1]
    range_start = rangefield[2]
    range_end = rangefield[3]


##############################################################################
# MAIN PROGRAM
##############################################################################

# Parse command line
if len(sys.argv) != 3:
    if len(sys.argv) == 1 && sys.argv[1] == 'help':
        sys.stdout.write(USAGE)
        sys.stdout.write(HELP)
        sys.exit(0)
    sys.stderr.write(USAGE)
    sys.exit(1)
original_sequence_filename = sys.argv[1]
parameter_str = sys.argv[2]

# Open Sequence
get_sequence(sys.argv[1])

# Get parameters
params = validate_parameters(parameter_str.split(','))

# Print sequence info
sys.stderr.write("Sequence header = %s \n" % seq_id)
sys.stderr.write("Sequence length = %d \n" % (len(SEQUENCE)))

# Generate sequences
SeqList = generate_sequences(Bases, Spacerlen, Spacerpos, Begpos, Endposition)

#resort list
SeqList = sortlist(SEQUENCE, SeqList)

# get rid of original sequence
del SeqList[:1]

# Number of sequences generated
sys.stderr.write("Number of sequences made = %d \n" % (len(SeqList)-1))

#get header info
header = seq_id.split()
gene_id = '>' + header[0]
rangefield = re.split(HEADER_SPLIT, header[1])
chrom_id = rangefield[1]
range_start = rangefield[2]
range_end = rangefield[3]

finalout = []

#get index of changes
for key, seq in enumerate(SeqList):
    changes = [i for i, (s1, s2) in enumerate(zip(SEQUENCE,seq)) if s1 != s2]
    finalout.append({"key":key+1, "mutations start": changes[0]+1 , "mutations":len(changes[0:]), "mutations end":changes[-1]+1, "sequence": seq})

#sort by index of change
finalout = sorted(finalout, key=itemgetter("mutations start", "mutations end", "mutations"))

finalout = {"gene id": gene_id, "chrom": chrom_id, "range start": range_start, "sequence length": len(SEQUENCE),\
             "base sequence": {"sequence": SEQUENCE}, "range end": range_end, "number of sequences": len(finalout), "list": finalout}
