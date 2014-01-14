BASES = ['A','C','G','T']
NUM_BASES = len(BASES)

def mutateseq(seq, pos):
    nseq = seq[:pos]
    nseq += BASES[(BASES.index(seq[pos])+1)%4]
    nseq += seq[pos+1:]
    return nseq

def mutatebase(base):
    nbase = BASES[(BASES.index(base)+1)%4]
    return nbase
    
def permutatebase(base):
    b1 = mutatebase(base)
    b2 = mutatebase(b1)
    b3 = mutatebase(b2)
    return [b1,b2,b3]

def permutateseq(seq, pos):
    seq1 = mutateseq(seq,pos)
    seq2 = mutateseq(seq1,pos)
    seq3 = mutateseq(seq2,pos)

    lseq = [[seq1], [seq2], [seq3]]
    return lseq
   
def permuteAll(seq):
    num_comb = NUM_BASES**len(seq)
    retlist = [""]*num_comb
    for curr_comb in range(num_comb):
        for curr_pos in range(len(seq)):
            retlist[curr_comb] += BASES[curr_comb/(NUM_BASES**curr_pos)%NUM_BASES]    
    retlist = sortlist(retlist, seq)
    
    return retlist

def generate(seq, bases=1, spacelen=0, spacepos=1, beg=1, end=-1):
    #set end to length of seq
    if end == -1:
        end = len(seq)

    #check for parameter errors
    if bases<1 or bases>end:
        print "bases parameter error"
        return
    if beg<1 or beg>end:
        print "beginning position error"
        return
    if end<beg or end>len(seq):
        print "end position error"
        return
    if bases>1:
        if spacelen<0 or spacelen>(bases-2):
            print "spacer length error"
            return
        if spacepos<=0 or spacepos+spacelen>=bases:
            print "spacer position error"
            return
    else:
        if spacelen!=0:
            print "spacer length error"
            return
        if spacepos!=1:
            print "spacer position error"
            return
    
    #fix positions to match array 0-n 
    beg -= 1
    
    #return list with permutations
    retlist = [seq]
    
    #initialize while loop variables
    #first slice initialize
    slicebeg = beg      #beginning of slice in sequence
    sliceend = slicebeg+bases   #end of slice in sequence
    spacerend = spacepos+spacelen
    
    #get permutations
    while (sliceend <= end):
        #get sequence parts before and after slice
        preslice = seq[:slicebeg]
        postslice = seq[sliceend:]
        permlist = []
        #get slice 
        seqslice = seq[slicebeg:sliceend]
        #part of slice that gets permutated
        permslice = seqslice[:spacepos] + seqslice[spacerend:]
        #spacer slice
        spacer = seqslice[spacepos:spacerend]
        #get permutations
        perms = permuteAll(permslice)
        #get rid of duplicate starting sequence in permutated list
        del perms[0]
        
        #reattach permutations with spacer
        for perm in perms:
            connectedslice = perm[:spacepos] + spacer + perm[spacepos:]
            permlist.append(connectedslice)
        #put slice back into sequence parts
        for perm in permlist:
            retlist.append(preslice + perm + postslice)    
        
        #increment slice start and end positions    
        slicebeg += 1
        sliceend = slicebeg+bases
    
    #eliminate duplicates    
    retlist = checkduplicates(retlist)
    #resort list
    sortlist(retlist, seq)
    
    return retlist
    
def checkduplicates(lis):   
    seen = set()
    result = []
    for item in lis:
        if item not in seen:
            seen.add(item)
            result.append(item)
    
    return result

def sortlist(list0, seq):
    ind = list0.index(seq)
    retlist = list0[ind:] + list0[:ind]
    return retlist
  
def print_nestedlist(lis):
    for item in lis:
        if isinstance(item, list):
            printlist(item)
        else:
            print item

def printlist(lis, cols=3):
    print lis[0], ' ' , len(lis)
    for i in range(1,len(lis)):
        print lis[i], ' ',
        if i%cols == 0:
            print ''

seq = "AAAAAAAAAAAA"

gene = generate(seq,4,2,1) 

printlist( gene)