import subprocess
import sys
import ConfigParser
import json
from operator import itemgetter


km = 3
sp = 0
CP = ConfigParser.ConfigParser()
CP.read('config.cfg')
lowestr_file = open("outputs/lowest.txt", "w+b")
js_lowlist = []
 
for i in range(18):
    ckm = km+i
    csp = sp+i
    print i, ckm, csp
    CP.set('Params', 'kmers', ckm)
    CP.set('Params', 'spacer length', csp)
    CP.set('Options', 'rvalues file', 'outputs/rvalues/' + str(ckm) + '-' + str(csp) + '_rvalues.txt')
    with open('config.cfg', 'wb') as configfile:
        CP.write(configfile)
    
    subprocess.call(['batch.py'], shell=True)
    subprocess.call(['graph.py'], shell=True)
    
    js = json.loads(open("outputs/select_data.txt", "r").read())
    js_lowlist.append(js["list"][0])
    
    
js_lowlist = sorted(js_lowlist, key=itemgetter("rval"))
lowestr_file.write(json.dumps(js_lowlist, indent=4))
lowestr_file.close()