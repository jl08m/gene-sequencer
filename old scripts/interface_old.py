from Tkinter import *
from ttk import Combobox, Frame, Progressbar
import tkFileDialog
import subprocess
import ConfigParser
import os
import time
MODEL_DIR = "./models/"
BATCH_FILE = ""

def openfile():
    input_file.set(tkFileDialog.askopenfilename())
    
    return

def opendir():
    output_dir.set(tkFileDialog.askdirectory())
    
    return
    
def setdefault():
    config = ConfigParser.ConfigParser()
    config.read('config.cfg')
    input_file.set(config.get('Files', 'input file'))
    output_dir.set(config.get('Files', 'output directory'))
    models_cb.set(config.get('Files', 'model file')[len(MODEL_DIR):])
    kmers.set(config.get('Params', 'kmers'))
    spacel.set(config.get('Params', 'spacer length'))
    spacep.set(config.get('Params', 'spacer after position'))
    startp.set(config.get('Params', 'starting at position'))
    endp.set(config.get('Params', 'ending at position'))
    cgraphs.set(config.get('Graph', 'create graphs'))
    topp.set(config.get('Graph', 'top'))
    botp.set(config.get('Graph', 'bottom'))
    kpgdata.set(config.getboolean('Options', 'keep graph data'))
    kprval.set(config.getboolean('Options', 'keep all rvalues'))
    kpseqs.set(config.getboolean('Options', 'keep all sequences'))
    
    return

def setconfig():
    config = ConfigParser.ConfigParser()
    config.read('config.cfg')
    
    config.set('Svm', 'input file', chrom_entry.get())
    config.set('File', 'output directory', output_entry.get() + '/')
    config.set('Svm', 'model file', MODEL_DIR + models_cb.get())
    
    config.set('List', 'kmers', kmers.get())
    config.set('List', 'spacer length', spacel.get())
    config.set('List', 'spacer after position', spacep.get())
    config.set('List', 'starting at position', startp.get())
    config.set('List', 'ending at position', endp.get())
    
    config.set('Graph', 'create graphs', cgraphs.get())
    config.set('Graph', 'top', topp.get())
    config.set('Graph', 'bottom', botp.get())
    
    config.set('List', 'keep all sequences', kpseqs.get())
    config.set('Rvals', 'keep rvals', kprval.get())
    config.set('Graph', 'keep graph data', kpgdata.get())
    
    with open('config.cfg', 'wb') as configfile:
        config.write(configfile)
    
    return 

def startall():
    setconfig()
    diag.insert(END, "Starting...\n")
    sp = subprocess.Popen(['batch.py'], shell=True, stdout=subprocess.PIPE)
    diag.insert(END, sp.communicate()[0])
    return

def stopall():
    return

root = Tk()
root.title("Main")
 

#VARIABLES------------------------------
input_file = StringVar()
output_dir = StringVar()

model = StringVar()

kmers = IntVar()
spacel = IntVar()
spacep = IntVar()
startp = IntVar()
endp = IntVar()

cgraphs = BooleanVar()
topp = IntVar()
botp = IntVar()

kpseqs = BooleanVar()
kprval = BooleanVar()
kpgdata = BooleanVar()

#CHILD FRAMES------------------------------
rootframe = Frame(root, padding="10 10 10 10")


mainframe = Frame(rootframe)
mainframe.grid(column=0, row=0, sticky=(N, W, E, S))
mainframe.columnconfigure(0, weight=1)
mainframe.rowconfigure(0, weight=1)

for child in mainframe.winfo_children(): child.grid_configure(padx=5, pady=5)

fileframe = Frame(mainframe)
fileframe.grid(column=1, row=1, sticky=(W, E, N), columnspan=4)

paramsframe = LabelFrame(mainframe, text="Parameters")
paramsframe.grid(column=1, row=3, sticky=(W, N, E))

graphframe = LabelFrame(mainframe, text="Graph options")
graphframe.grid(column=2, row=3, sticky=(W, N, E))

optionsframe = LabelFrame(mainframe, text="Output Options")
optionsframe.grid(column=3, row=3, sticky=(W, N, E))

dialogframe = Frame(rootframe)
statusframe = Frame(rootframe)


#FILE SEARCH BLOC------------------------------

chrom_entry = Entry(fileframe, width=30, textvariable=input_file)
output_entry = Entry(fileframe, width=30, textvariable=output_dir)
chrom_button = Button(fileframe, text="Find file", command=openfile)
output_button = Button(fileframe, text="Find directory", command=opendir)

chrom_button.grid(column=1, row=1, sticky=(W, E))
chrom_entry.grid(column=2, row=1, sticky=(W, E), columnspan=4)
output_button.grid(column=1, row=2, sticky=(W, E))
output_entry.grid(column=2, row=2, sticky=(W, E), columnspan=4)

#MODEL BLOCK------------------------------
models_lb = Label(fileframe, text="Models")
models_cb = Combobox(fileframe, width=20, state='readonly')
modellist = []
for files in os.listdir(MODEL_DIR):
    if files.endswith(".txt"):
        modellist.append(files)
models_cb['values'] = modellist
models_lb.grid(column=1, row=3, sticky=(W))
models_cb.grid(column=2, row=3, sticky=(W))


#PARAMETER BLOCK------------------------------
kmer_label = Label(paramsframe, text="Kmers")
spacerl_label = Label(paramsframe, text="Spacer Length")
spacerp_label = Label(paramsframe, text="Spacer after position")
startpos_label = Label(paramsframe, text="Starting at position")
endpos_label = Label(paramsframe, text="Ending at position")

kmer_entry = Entry(paramsframe, width=5, textvariable=kmers)
spacerl_entry = Entry(paramsframe, width=5, textvariable=spacel)
spacerp_entry = Entry(paramsframe, width=5, textvariable=spacep)
startpos_entry = Entry(paramsframe, width=5, textvariable=startp)
endpos_entry = Entry(paramsframe, width=5, textvariable=endp)

kmer_label.grid(column=1, row=1, sticky=(W))
spacerl_label.grid(column=1, row=2, sticky=(W))
spacerp_label.grid(column=1, row=3, sticky=(W))
startpos_label.grid(column=1, row=4, sticky=(W))
endpos_label.grid(column=1, row=5, sticky=(W))

kmer_entry.grid(column=2, row=1, sticky=(W,E))
spacerl_entry.grid(column=2, row=2, sticky=(W,E))
spacerp_entry.grid(column=2, row=3, sticky=(W,E))
startpos_entry.grid(column=2, row=4, sticky=(W,E))
endpos_entry.grid(column=2, row=5, sticky=(W,E))


#GRAPH BLOCK------------------------------
graph_check = Checkbutton(graphframe, text="Create Graphs", variable=cgraphs, onvalue=True, offvalue=False)
top_lb = Label(graphframe, text="Top %")
bot_lb = Label(graphframe, text="Bottom %")
top_entry = Entry(graphframe, width=5, textvariable=topp)
bot_entry = Entry(graphframe, width=5, textvariable=botp)

graph_check.grid(column=0, row=1, sticky=(W), columnspan=2)
top_lb.grid(column=0, row=2, sticky=(W))
bot_lb.grid(column=0, row=3, sticky=(W))
top_entry.grid(column=1, row=2, sticky=(W))
bot_entry.grid(column=1, row=3, sticky=(W))


#OPTIONS BLOCK------------------------------
keepseqs_check = Checkbutton(optionsframe, text="Keep all sequences", variable=kpseqs, onvalue=True, offvalue=False)
keeprvals_check = Checkbutton(optionsframe, text="Keep all r-values", variable=kprval, onvalue=True, offvalue=False)
keepgdata_check = Checkbutton(optionsframe, text="Keep all graphed data", variable=kpgdata, onvalue=True, offvalue=False)

keepseqs_check.grid(column=0, row=1, sticky=(W))
keeprvals_check.grid(column=0, row=2, sticky=(W))
keepgdata_check.grid(column=0, row=3, sticky=(W))

#------------------------------
defaults_bt = Button(mainframe, text="Last Settings", command=setdefault)
defaults_bt.grid(column=3, row=1)

start_bt = Button(statusframe, text="Start", width=10, command=startall)
start_bt.pack(side=LEFT)
progress = Progressbar(statusframe, orient=HORIZONTAL, length=450, mode='determinate')
progress.pack(side=LEFT)
stop_bt = Button(statusframe, text="Stop", width=10, command=stopall)
stop_bt.pack(side=LEFT)

diag = Text(dialogframe)
scroll = Scrollbar(dialogframe)
scroll.pack(side=RIGHT, fill=Y)
diag.pack(side=LEFT, fill=Y)
scroll.config(command=diag.yview)
diag.config(yscrollcommand=scroll.set)
#------------------------------
rootframe.pack(fill=BOTH)
mainframe.pack()
dialogframe.pack()
statusframe.pack()
setdefault()

#diag.insert(END, sys.stdout)

root.mainloop()