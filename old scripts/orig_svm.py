import subprocess
import ConfigParser

cf = ConfigParser.ConfigParser()
cf.read('config.cfg')

svm_chrom = cf.get('Resources', 'svm file')
input = cf.get('Files', 'input file')
model = cf.get('Files', 'model file')
svm_file = cf.get('Orig', 'orig svms')

svms = open(svm_file, "w+b")

subprocess.call([svm_chrom, input, model], shell=True, stdout=svms)