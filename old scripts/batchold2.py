import sys
import subprocess
import tempfile
from os import unlink
from time import time
from operator import itemgetter
from scipy.stats.stats import pearsonr

##argsfile = sys.argv[1]
#argsfile = "args.txt"
#
#args = []
#for line in open(argsfile, "r+b").readlines():
#    args.append(line.strip())
#
#cmd_svm = args[0]
#cmd_gen = args[1]
#chrom_filename = args[2]
#model_filename = args[3]
#params = args[4]
#seqlist_filename = args[5]
#rvals_filename = args[6]
#msgs_filename = args[7]

svm_filename = "svm-chrom.py"
gen_filename  = "genseq.py"
input_chrom_filename  = "fasta/test_atm.fasta"
model_filename  = "a375_nuc.model.txt"
full_seqlist_filename  = "full_seqlist.txt"
rvalues_filename  = "rvals.txt"
parameters = "5,3,1,60,70"
#Open files for writing
full_seqlist_file = open(full_seqlist_filename, "w+b")
rvalues_file = open(rvalues_filename, "w+b")
full_seqlist_file = open(full_seqlist_filename, "w+b")

# Get SVM values for original sequence
temp_orig_svm_file = tempfile.TemporaryFile()
subprocess.call([svm_filename, input_chrom_filename, model_filename], shell=True, stdout=temp_orig_svm_file, stderr=tempfile.TemporaryFile())
temp_orig_svm_file.seek(0)
orig_svmvals = []
# Retrieve svm values
for line in temp_orig_svm_file.readlines():
    orig_svmvals.append(float(line.strip()))
    

# Generate list of Sequences
subprocess.call([gen_filename, input_chrom_filename, parameters], shell=True, stdout=full_seqlist_file, stderr=tempfile.TemporaryFile())
full_seqlist_file.seek(0)
full_seqlist = full_seqlist_file.readlines()

starttime = time()

rval_list = []
# For each generated sequence, create temp file that has name and sequence fasta style
for i,seq in enumerate(full_seqlist):
    #separate name from sequence
    p = seq.split('|')
    seq_id = p[0]
    seq_seq = p[1]
    
    # Write sequence into temp file 
    temp_seq_file = tempfile.NamedTemporaryFile(delete=False)
    temp_seq_file.write(seq_id + '\n')
    temp_seq_file.write(seq_seq)

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
        
#    #test print svm vals 
#    print 'origsvm= ', orig_svmvals
#    print 'tempsvm= ', temp_svmvals
    
    # Compare SVMs and get r-Values
    rval = pearsonr(orig_svmvals, temp_svmvals)
    
    # put rval with seq_id into output file
#    rvalues_file.write(seq_id + ' | ' + str(rval) + '\n')
    rval_list.append((seq_id, rval))
    
    # unlink temp file for deletion
    temp_seq_file.close()
    unlink(temp_seq_file.name)

    sys.stderr.write("\r%d / %d" %(i, len(full_seqlist)))
    sys.stderr.flush()

#    stop after #
#    if i==0:
#        break

endtime = time()
sys.stderr.write("\ntime to complete: %d seconds" %(endtime-starttime))
#sys.stderr.write("\nfor %d sequences" %i)
rval_list = sorted(rval_list, key=itemgetter(1), reverse=True)

for item in rval_list:
    r = str(item[1])
    ent = item[0] + ' | ' + r
    rvalues_file.write(ent + '\n')
sys.stderr.write("\ntime to complete: %d seconds" %(endtime-starttime))