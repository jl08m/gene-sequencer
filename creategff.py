import json
import getconfigs
import sys

def converttogff(chrom, position, svms):
    retlist = []
    V1 = chrom
    V2 = '.'
    V3 = '.'
    V4 = '.' #start position
    V5 = '.' #end position
    V6 = '.' #svm
    V7 = '.'
    V8 = '.'
    V9 = ";color=000000\n"
    t = '\t'
    v1 = V1 + t + V2 + t + V3 + t
    v2 = V7 + t + V8 + t + V9
    for i,svm in enumerate(svms):
        vpos = str(int(position)+i+26)
        retlist.append(v1 + vpos + t + vpos + t + str(svm) + t + v2)

    return retlist

###############################################

GRAPH_DIR = "outputs/graph_data/"
GFF_DIR = "gffs/"
cf = getconfigs.configs()

#get relevant data
selected_data = open(cf.selectedfile)
json_all = json.loads(selected_data.read())
json_baseseq = json_all["base sequence"]
json_list = json_all["list"]
chrom = json_all["chrom"]
rangestart = json_all["range start"]

##base sequence gff part
#list of gff lines for base sequence
baseseq_gff = converttogff(chrom, rangestart, json_baseseq["svmvals"])
#open base gff file for writing
baseseq_gff_file = open(GRAPH_DIR+"baseseq.gff", "w+b")
#write into gff file
for line in baseseq_gff:
    baseseq_gff_file.write(line)
#close file
baseseq_gff_file.close()

for i, seq in enumerate(json_list): 
    current_seq_svms= seq["svmvals"]
    cmpgff = converttogff(chrom, rangestart, current_seq_svms)
    cmp_filename = GRAPH_DIR + GFF_DIR + \
        str(i) + '_' + json_all["gene id"][1:] + '_' + \
        str(seq["mutations start"]) + '-' + str(seq["mutations end"]) +\
        '_' + str(seq["key"]) + '.gff'
        
#    print cmp_filename
    cmpgff_gff_file = open(cmp_filename, "w+b") 
#    tempfile = tempfile.NamedTemporaryFile(delete=False)
#    for line in cmplist:
#        tempfile.write(line)
#        print line,
#    
    #create gff plot
    # do R stuff
    for line in cmpgff:
        cmpgff_gff_file.write(line)
    
#    tempfile.close()
#    unlink(tempfile.name)
#    break
    cmpgff_gff_file.close()
    
sys.stderr.write("created gffs\n")
