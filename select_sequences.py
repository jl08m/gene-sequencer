#Written by Jia Li
#Florida State University, Department of Biology
#07-25-2012

#################################
#This script takes the list of sequences sorted by r-value
#then generates a new list of the top and bottom percentages
#specified in the configuration file
#################################

import tempfile
import sys
import json
from os import unlink
import getconfigs



#get config stuff
cf = getconfigs.configs()

rval_filename = cf.rvalfile
selectlist_filename= cf.selectedfile
gettop = cf.topcut
getbot = cf.botcut

#open files
rval_file = open(rval_filename,"r+b")
selectedlist_file = open(selectlist_filename, "w+b")

#set original svm value array
rval_list = json.loads(rval_file.read())
svm_list = []
for line in rval_list["base sequence"]["sequence"]:
    svm_list.append(line)

#intialize json output
json_out = {}
json_out["top %"] = str(gettop)
json_out["bottom %"] = str(getbot)
json_out["gene id"] = rval_list["gene id"]
json_out["base sequence"] = rval_list["base sequence"]
json_out["chrom"] = rval_list["chrom"]
json_out["range start"] = rval_list["range start"]

#get top and bottom cutoff sequences list  
listsize = len(rval_list["list"]) 
top = []
bot = []

#get top percentage
if gettop >= 0:
    topper = float(gettop) / 100.0
    topcut = int(round(topper*listsize))
    top = rval_list["list"][:topcut]
#    json_out["top"] = top
    
#get bottom percentage
if getbot >= 0:
    botper = float(getbot) / 100.0
    botcut = int(round(botper*listsize))
    bot = rval_list["list"][-botcut:]
#    bot.reverse()
#    json_out["bottom"] = bot
   
#combine into 1 list
json_out["list"] = top + bot
    
#message to show how many were gotten
sys.stderr.write("top %d%%: %d sequences\n" %(gettop,len(top)))
sys.stderr.write("bottom%d%%: %d sequences\n" %(getbot,len(bot)))

#write to output file
selectedlist_file.write(json.dumps(json_out, indent=2))
t = open('temp/top.txt', 'w')
t.write(json.dumps(top, indent=2))
b = open('temp/bot.txt', 'w')
b.write(json.dumps(bot, indent=2))


