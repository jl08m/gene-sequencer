#Written by Jia Li
#Florida State University, Department of Biology
#07-25-2012

#################################
#Creates a file of all r-values with the
#relevant information for each sequence generated
#There are a lot of code commented out
#used for testing purposes
#################################

import sys
import subprocess
import tempfile
from os import unlink
from time import time
from operator import itemgetter
from scipy.stats.stats import pearsonr
import json
import getconfigs

#the output list scrambles the order of the list of dicts for some reason
#use Ordereddicts if the output text file is important


cf = getconfigs.configs()

svm_filename = cf.svm_script
gen_filename = cf.gen_script
input_chrom_filename = cf.input_file
model_filename = cf.model
full_seqlist_filename = cf.seqfile
rvalues_filename = cf.rvalfile

#set parameters
parameters = str(cf.kmers) + ',' + str(cf.spacelen) + ',' + \
    str(cf.spacepos) + ',' + str(cf.startpos) + ',' + str(cf.endpos)

#notification messages 
sys.stderr.write(full_seqlist_filename + '\n')
sys.stderr.write(rvalues_filename + '\n')

#Open files for writing
rvalues_file = open(rvalues_filename, "w+b")

# Get SVM values for original sequence
orig_svmvals=[]
###################################
#if no original svm values are saved then
#use this section to generate a temporary file
#for the use in this script
###################################
temp_orig_svm_file = tempfile.TemporaryFile()
subprocess.call([svm_filename, input_chrom_filename, model_filename], shell=True, stdout=temp_orig_svm_file, stderr=sys.stdout)
temp_orig_svm_file.seek(0)
# Retrieve svm values
for line in temp_orig_svm_file.readlines():
    orig_svmvals.append(float(line.strip()))
#otherwise use these lines
#for line in orig_svms.readlines():
#    orig_svmvals.append(float(line.strip()))
    

# Generate list of Sequences
subprocess.call([gen_filename, input_chrom_filename, parameters, full_seqlist_filename], shell=True, stderr=sys.stdout)
json_seqlist = json.loads(open(full_seqlist_filename, "r+b").read())
#print json.dumps(json_seqlist, indent=2)
starttime = time()

#add base sequence to full rvalues json
json_seqlist["base sequence"]["svmvals"] = orig_svmvals

#initialize lists
rval_list = []
listsize = len(json_seqlist["list"])

# For each generated sequence, create temp file that has name and sequence fasta style
for i, seq in enumerate(json_seqlist["list"]):
    
    #CREATE TEMPORARY FILE FOR THIS SEQUENCE
    #get name for header
    seq_head = json_seqlist["gene id"] + \
                    str(seq["key"]) + \
                    str(seq["mutations start"]) + \
                    str(seq["mutations"]) + \
                    str(seq["mutations end"])
                    
    # Write sequence into temp file 
    temp_seq_file = tempfile.NamedTemporaryFile(delete=False)
    temp_seq_file.write(seq_head + '\n')
    temp_seq_file.write(seq["sequence"] + ' ')

## test seq file, print content
#    temp_seq_file.seek(0)
#    for line in temp_seq_file.readlines():
#        print line,

    # Get svm values and put into temp file
    temp_seq_svmvals_file = tempfile.TemporaryFile()
    temp_seq_file.seek(0)
    subprocess.call([svm_filename, temp_seq_file.name, model_filename], shell=True, stdout=temp_seq_svmvals_file, stderr=tempfile.TemporaryFile())
    temp_seq_svmvals_file.seek(0)
    
#  # test seq svm file, print content
#    temp_seq_svmvals_file.seek(0)
#    print 'tempseqsvmfile'
#    for line in temp_seq_svmvals_file.readlines():
#        print line,
  
    # Retrieve svm values
    temp_svmvals = []
    temp_seq_svmvals_file.seek(0)
    for line in temp_seq_svmvals_file.readlines():
        temp_svmvals.append(float(line.strip()))
        
    #put svm values into the final JSON also
    #WARNING CREATES LONG LIST OF NUMBERS
    seq["svmvals"] = temp_svmvals
    
#    #test print svm vals 
#    print 'origsvm= ', orig_svmvals
#    print 'tempsvm= ', temp_svmvals
    
    # Compare SVMs and get r-Values
    rval = pearsonr(orig_svmvals, temp_svmvals)
    
    # put rval with seq_id into output file
    #rval_list.append((seq_id, rval))
#    print 'rval=', rval
    seq["rval"] = rval
    
    # unlink temp file for deletion
    temp_seq_file.close()
    unlink(temp_seq_file.name)

    progress = int(round(float(i)/float(listsize)*float(100)))
    #sys.stderr.write("%d %%\n" %progress)
    sys.stdout.write('%d/%d\n' % (float(i), float(listsize)))

#    stop after #
#    if i==0:
#        break

endtime = time()
#sys.stderr.write("\nfor %d sequences" %i)
#print 'jsonseqlist=' ,json_seqlist["list"][0]
#json_seqlist["number of sequences"] = len(json_seqlist["list"])
json_seqlist["list"] = sorted(json_seqlist["list"], key=itemgetter("rval"), reverse=True)

#write and close the rvalue file
rvalues_file.write(json.dumps(json_seqlist, indent=2))
rvalues_file.close()

#time to completion
sys.stdout.write("\ntime to complete: %d seconds\n" % (endtime - starttime))
