import itertools

class Sequence:
    
    #if more bases used, add to list
    BASES = ['A','C','G','T']
    NUM_BASES = len(BASES)
    
        
##############################################################################
    #constructor for class
    def __init__(self, arg1=None, arg2=None):
        #no argument defaults to length 10 all 'A'
        if arg1 is None:
            self.sequence = "AAAAAAAAAA"
        #if int is parameter than make int length
        elif isinstance(arg1, int):
            #if 2nd argument, then make length of that base
            if arg2 == 'C' or arg2 == 'G' or arg2 == 'T': 
                self.sequence = arg2 * arg1
            #else make length of A    
            else:
                self.sequence = 'A' * arg1
        #if string then make sequence that string 
        elif isinstance(arg1, str):
            self.sequence = arg1
            
        else:
            self.sequence = None
            print "Error"
        
        self.seqlist = [self.sequence]
    
##############################################################################
    #print original sequence
    def printseq(self):
        print self.sequence
        
##############################################################################
    #print list in number of columns
    def printlist(self, cols=4):
        for i in range(len(self.seqlist)):
            if i%cols == 0:
                print ''
            print self.seqlist[i], ' ',
        print '\n'
            
##############################################################################
    #deprecated permutation function, use if itertools unavailable
    def permuteAll(self, seq):
        num_comb = self.NUM_BASES**len(seq)
        retlist = [""]*num_comb
        for curr_comb in range(num_comb):
            for curr_pos in range(len(seq)):
                retlist[curr_comb] += self.BASES[curr_comb/(self.NUM_BASES**curr_pos)%self.NUM_BASES]
    
        return retlist

##############################################################################
    #generate permutations and create list of sequences
    def generate(self, bases=1, spacelen=0, spacepos=1, beg=1, end=-1):
        #set end to length of sequence
        if end == -1:
            end = len(self.sequence)

        #check for parameter errors
        if not self.errorcheck(bases,spacelen,spacepos,beg,end):
            return 
        
        #fix positions to match array 0-n 
        beg -= 1
        
        #initialize while loop variables
        #first slice initialize
        slicebeg = beg      #beginning of slice in sequence
        sliceend = slicebeg+bases   #end of slice in sequence
        spacerend = spacepos+spacelen
        
        #get permutations
        while (sliceend <= end):
            #get sequence parts before and after slice
            preslice = self.sequence[:slicebeg]
            postslice = self.sequence[sliceend:]
            permslicelist = []
            
            #get slice 
            seqslice = self.sequence[slicebeg:sliceend]
            
            #part of slice that gets permutated
            permslice = seqslice[:spacepos] + seqslice[spacerend:]
            
            #spacer slice
            spacer = seqslice[spacepos:spacerend]
            
            #get permutations
            temp = itertools.product(self.BASES, repeat=len(permslice))
            perms = []
            for i in temp:
                perms.append("".join(i))
            
            #deprecated permutation part
            ##perms = self.permuteAll(permslice)
            ##get rid of duplicate starting sequence in permutated list
            ##del perms[0]
            
            #reattach permutations with spacer
            for perm in perms:
                permslicelist.append(perm[:spacepos] + spacer + perm[spacepos:])
                
            #put slice back into sequence parts
            for perm in permslicelist:
                self.seqlist.append(preslice + perm + postslice)    
            
            #increment slice start and end positions    
            slicebeg += 1
            sliceend = slicebeg+bases
        
        #eliminate duplicates    
        self.seqlist = self.checkduplicates()
        
        #resort list
        self.seqlist = self.sortlist()
    
##############################################################################
    #error parameters
    def errorcheck(self, bases, spacelen, spacepos, beg, end):
        ErrorCode = -1
        if bases<1 or bases>end:
            ErrorCode = 1
        elif beg<1 or beg>end:
            ErrorCode = 2
        elif end<beg or end>len(self.sequence):
            ErrorCode = 3
        elif bases>1:
            if spacelen<0 or spacelen>(bases-2):
                ErrorCode = 4 
            if spacepos<=0 or spacepos+spacelen>=bases:
                ErrorCode = 5
        else:
            if spacelen!=0:
                ErrorCode = 4 
            if spacepos!=1:
                ErrorCode = 5
                
        #testing line
        #ErrorCode = -1        
        
        if ErrorCode == 0:
            print "Unknown error"    
        elif ErrorCode == 1:
            print "bases parameter error"
        elif ErrorCode == 2:
            print "beginning position error"
        elif ErrorCode == 3:
            print "end position error"
        elif ErrorCode == 4:
            print "spacer length error"
        elif ErrorCode == 5:
            print "spacer position error"
            

        if ErrorCode == -1:
            return True 
        else:
            return False 
    
##############################################################################
    #generate all permutations for entire sequence
    def generateall(self):
        self.generate(len(self.sequence))
        
##############################################################################
    #check duplicates in list, returns list
    def checkduplicates(self):   
        seen = set()
        retlist = []
        for item in self.seqlist:
            if item not in seen:
                seen.add(item)
                retlist.append(item)
        return retlist         
                
##############################################################################
    #sort list by starting sequence
    def sortlist(self):
        ind = self.seqlist.index(self.sequence)
        retlist = self.seqlist[ind:] + self.seqlist[:ind]
        return retlist
    
##############################################################################
    #return size of list
    def size(self):
        return len(self.seqlist)
    
##############################################################################
    #return length of original sequence
    def length(self):
        return len(self.sequence)
    
##############################################################################
    #delete list of sequences
    def clearlist(self):
        del self.seqlist[1:]
        
##############################################################################