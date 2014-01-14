import subprocess
import getconfigs
import sys


cf = getconfigs.configs()
    
subprocess.call([cf.batch_script], shell=True, stdout=sys.stdout, stderr=sys.stderr)
subprocess.call([cf.select_script], shell=True, stdout=sys.stdout, stderr=sys.stderr)

if cf.creategraphs == True:
    subprocess.call([cf.creategff_script], shell=True, stdout=sys.stdout, stderr=sys.stderr)

sys.exit()