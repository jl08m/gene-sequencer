# FILE: genseq.py
# AUTHOR: Jia Li
# CREATE DATE: 2012-07-10
# VERSION: 1.0

import itertools
import sys

#if more bases used, add to list
BASES = ['A', 'C', 'G', 'T']
NUM_BASES = len(BASES)
Sequence = ""
SeqList = []

# example usage: geneseq.py somesequence.fasta 4,2,2,0,0
usage = """USAGE: genseq.py <sequence> <bases,spacer length, spacer position, beginning, end>\n"""
    
#############################################################################
#generate permutations and create list of sequences
def generate(bases=1, spacelen=0, spacepos=1, beg=1, end=-1):
    retlist = []
    #set end to length of sequence
    if end == -1:
        end = len(Sequence)
    
    #fix positions to match array 0-n 
    beg -= 1
    
    #initialize while loop variables
    #first slice initialize
    slicebeg = beg      #beginning of slice in sequence
    sliceend = slicebeg + bases   #end of slice in sequence
    spacerend = spacepos + spacelen
    
    #get permutations
    while (sliceend <= end):
        #get sequence parts before and after slice
        preslice = Sequence[:slicebeg]
        postslice = Sequence[sliceend:]
        permslicelist = []
        
        #get slice 
        seqslice = Sequence[slicebeg:sliceend]
        
        #part of slice that gets permutated
        permslice = seqslice[:spacepos] + seqslice[spacerend:]
        
        #spacer slice
        spacer = seqslice[spacepos:spacerend]
        
        #get permutations
        temp = itertools.product(BASES, repeat=len(permslice))
        perms = []
        for i in temp:
            perms.append("".join(i))
        
        #deprecated permutation part
        ##perms = permuteAll(permslice)
        ##get rid of duplicate starting sequence in permutated list
        ##del perms[0]
        
        #reattach permutations with spacer
        for perm in perms:
            permslicelist.append(perm[:spacepos] + spacer + perm[spacepos:])
            
        #put slice back into sequence parts
        for perm in permslicelist:
            retlist.append(preslice + perm + postslice)    
        
        #increment slice start and end positions    
        slicebeg += 1
        sliceend = slicebeg + bases
    

    return retlist

##############################################################################
#error parameters
def errorcheck(bases, spacelen, spacepos, beg, end):
    ErrorCode = 0

    
    if bases < 1 or bases > end:
        ErrorCode = 1
    elif beg < 1 or beg > end:
        ErrorCode = 2
    elif end < beg or end > len(Sequence):
        ErrorCode = 3
    elif bases > 1:
        if spacelen < 0 or spacelen > (bases - 2):
            ErrorCode = 4 
        if spacepos <= 0 or spacepos + spacelen >= bases:
            ErrorCode = 5
    else:
        if spacelen != 0:
            ErrorCode = 4 
        if spacepos != 1:
            ErrorCode = 5
            
    return ErrorCode 
   
##############################################################################
#check duplicates in list, returns list
def checkduplicates(seqlist):   
    seen = set()
    retlist = []
    for item in seqlist:
        if item not in seen:
            seen.add(item)
            retlist.append(item)
    return retlist         
            
##############################################################################
#sort list by starting sequence
def sortlist(sequence, seqlist):
    ind = seqlist.index(sequence)
    retlist = seqlist[ind:] + seqlist[:ind]
    return retlist

##############################################################################
# MAIN PROGRAM
##############################################################################

# Parse command line
if (len(sys.argv) != 3):
    sys.stderr.write(usage)
    sys.exit(1)
seq_filename = sys.argv[1]
parameters_str = sys.argv[2]


# Get parameters
parameters = parameters_str.split(',')
if len(parameters) != 5:
    sys.stderr.write("Incorrect Parameter\n")
    sys.exit(1)
Bases = int(parameters[0])
Spacerlen = int(parameters[1])
Spacerpos = int(parameters[2])
Begpos = int(parameters[3])
Endposition = int(parameters[4])

# Open the sequence file.
if (seq_filename == "-"):
    seq_file = sys.stdin
else:
    seq_file = open(seq_filename, "r")

# Read the first line.
seq_id = seq_file.readline()
seq_id = seq_id[1:-1]

# Read the sequence into memory
for line in seq_file:
    line = line[:-1]
    Sequence = Sequence + line

# Print sequence info
sys.stderr.write("Sequence length = %d \n" % (len(Sequence)))
sys.stderr.write("Sequence header = %s \n" % seq_id)

# check for parameter errors
if Bases == 0:
    Bases = 1
if Spacerlen == 0:
    Spacerlen = 0
if Spacerpos == 0:
    Spacerpos = 1
if Begpos == 0:
    Begpos = 1
if Endposition == 0:
    Endposition = len(Sequence)
    
error = errorcheck(Bases, Spacerlen, Spacerpos, Begpos, Endposition)
sys.stderr.write("Bases = %d \n" % (Bases))
sys.stderr.write("Spacer Length = %d \n" % (Spacerlen))
sys.stderr.write("Spacer Placed After Position %d \n" % (Spacerpos))
sys.stderr.write("Mutations started at position %d \n" % (Begpos))
sys.stderr.write("Mutations end at position %d \n" % (Endposition))
if error != 0:
    if error == 1:
        sys.stderr.write("bases parameter error \n")
    elif error == 2:
        sys.stderr.write("beginning position error \n")
    elif error == 3:
        sys.stderr.write("end position error \n")
    elif error == 4:
        sys.stderr.write("spacer length error \n")
    elif error == 5:
        sys.stderr.write("spacer position error \n")
    else:
        sys.stderr.write("Unknown error \n")
    sys.exit(1)
    
# Generate sequences
SeqList = generate(Bases, Spacerlen, Spacerpos, Begpos, Endposition)

#eliminate duplicates    
SeqList = checkduplicates(SeqList)

#resort list
SeqList = sortlist(Sequence, SeqList)

# get rid of original sequence
del SeqList[:1]

# Number of sequences generated
sys.stderr.write("Number of sequences made = %d \n" % (len(SeqList)-1))

# get indexes of differences for naming
diff = []
for q,seq in enumerate(SeqList):
    zz = [i for i, (s1, s2) in enumerate(zip(Sequence,seq)) if s1 != s2]
    diff.append(str(q+1) + ':' + str(zz[0]+1) + '-' + str(zz[-1]+1))

if len(SeqList) == 1:
    sys.stderr.write("Something wrong \n")
else:
    for i,j in zip(diff, SeqList):
        sys.stdout.write("%s:%s\n" %(i,j))
