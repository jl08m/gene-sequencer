
import re 

s = ">hg19_dna range=chr11:108093980-108094130 5'pad=0 3'pad=0 strand=+ repeatMasking=none"

t = re.split(' |=|:|-| ', s)
print t



