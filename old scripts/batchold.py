import sys
import subprocess
import tempfile
from os import unlink
from time import time
from scipy.stats.stats import pearsonr


#argsfile = sys.argv[1]
argsfile = "args.txt"

args = []
for line in open(argsfile, "r+b").readlines():
    args.append(line.strip())

cmd_svm = args[0]
cmd_gen = args[1]
chrom_filename = args[2]
model_filename = args[3]
params = args[4]
seqlist_filename = args[5]
rvals_filename = args[6]
msgs_filename = args[7]

msgs_file = open(msgs_filename, "w+b")
rvals_file = open(rvals_filename, "w+b")

# Get SVM values for original sequence
svmvalsFile = tempfile.TemporaryFile()
subprocess.call([cmd_svm, chrom_filename, model_filename], shell=True, stdout=svmvalsFile, stderr=msgs_file)
svmvalsFile.seek(0)
orig_svmvals = []
# Retrieve svm values
for line in svmvalsFile.readlines():
    orig_svmvals.append(float(line.strip()))
    

# Generate list of Sequences
seqlist_file = open(seqlist_filename,"w+b")
subprocess.call([cmd_gen, chrom_filename, params], shell = True, stdout=seqlist_file, stderr=msgs_file)
seqlist_file.seek(0)
seqlist = seqlist_file.readlines()

starttime = time()

# For each generated sequence, create temp file that has name and sequence fasta style
for i,seq in enumerate(seqlist):
    #separate name from sequence
    p = seq.split(':')
    seq_id = p[0]
    seq_seq = p[1]
    
    # Write sequence into temp file 
    tempfile_chrom = tempfile.NamedTemporaryFile(delete=False)
    tempfile_chrom.write('>' + seq_id + '\n')
    tempfile_chrom.write(seq_seq + ' ')

# test file, print content
#    tempfile_chrom.seek(0)
#    for line in tempfile_chrom.readlines():
#        print line,

    # Get svm values and put into temp file
    tempfile_svmvals = tempfile.TemporaryFile()
    tempfile_chrom.seek(0)
    subprocess.call([cmd_svm, tempfile_chrom.name, model_filename], shell=True, stdout=tempfile_svmvals, stderr=tempfile.TemporaryFile())
    tempfile_svmvals.seek(0)
    
    # Retrieve svm values
    temp_svmvals = []
    for line in tempfile_svmvals.readlines():
        temp_svmvals.append(float(line.strip()))
    
#    print orig_svmvals
#    print temp_svmvals
    
    # Compare SVMs and get r-Values
    rval = pearsonr(orig_svmvals, temp_svmvals)
    
    # put rval with seq_id into output file
    rvals_file.write(seq_id + ':' + str(rval) + '\n')
    
    tempfile_chrom.close()
    unlink(tempfile_chrom.name)

    sys.stderr.write("\r%d / %d" %(i, len(seqlist)))
    sys.stderr.flush()

#    stop after #
    if i==1:
        break

endtime = time()
sys.stderr.write("\ntime to complete: %d seconds" %(endtime-starttime))
sys.stderr.write("\nfor %d sequences" %i)
