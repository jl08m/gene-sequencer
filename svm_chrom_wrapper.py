__author__ = 'Jia'

import subprocess
import logging
CHROM_FILE = "svm-chrom.py"

logging.basicConfig(filename='log',
                    filemode='a',
                    format='%(asctime)s,%(msecs)d %(filename)s %(levelname)s %(funcName)s \n%(message)s',
                    datefmt='%H:%M:%S')


def svm_chrom(chromosome_file=None, model_file=None):
    """returns a list from the output of svm_chom.py.
    returns None if exception caught.

    Keyword arguments:
    chromosome_file -- the file containing the chromosome in the format:
        <{chromosome description}
        ACGTCGTA...
    model_file -- the file containing the model to compare against.
    """
    try:
        raw_string = subprocess.check_output([CHROM_FILE, chromosome_file, model_file], shell=True, stderr=open('tmp.log', 'a'))
        out = raw_string.split('\n')
        out.remove('')
        logging.debug(out)
        return out
    except subprocess.CalledProcessError, e:
        logging.error(e)
        return None