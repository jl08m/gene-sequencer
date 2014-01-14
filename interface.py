#Written by Jia Li
#Florida State University, Department of Biology
#07-25-2012

##################################
#This is the interface file for the suite of scripts
#Right now it should only allow the setting of the
#configurations file. Hopefully I will be able to figure
#out how to adapt some sort of progress bar and add more
#options.
##################################


from Tkinter import *
import ttk
import sys
import subprocess
import tkFileDialog
from time import sleep
from threading import Thread
from multiprocessing import Process
import os
import ConfigParser

class MainWindow:
    SCRIPT_FILE = 'script.py'
    CONFIG_FILE = 'config.cfg'
    MODEL_DIR = "models/"
   
   #initialize the class, set up all GUI widgets and variables 
    def __init__(self,parent):
       #widget variables
        self.inputsequence_path = StringVar()
        self.outputdirectory_path = StringVar()
        self.inputsequence_file = StringVar()
        self.outputdirectory_dir = StringVar()
            
        self.model = StringVar()
        
        self.kmers = IntVar()
        self.spacelen = IntVar()
        self.spacepos = IntVar()
        self.startpos = IntVar()
        self.endpos = IntVar()
        
        self.creategraphs = BooleanVar()
        self.topcut = IntVar()
        self.botcut = IntVar()
        
        self.keepseqs = BooleanVar()
        self.keeprvals = BooleanVar()
        
        self.progress = StringVar()
        
        
        #thread for batch
        self.script_thread = Thread(target=self.run_script)
#        self.batch_process = Process(target=self.run_batch)
        
        #initialize gui
        self.initfileframe(parent)
        self.initparamframe(parent)
        self.initselectframe(parent)
        self.initoptionsframe(parent)
        self.initmsgframe(parent)
        self.initbottombar(parent)  
         
        self.progress.set("0 %")
       
        #file selection frame
    def initfileframe(self, parent):
        self.fileframe = Frame(parent)
        self.fileframe.grid(column=0, row=0, sticky=W, columnspan=4)
        
        self.btn_opensequence = Button(self.fileframe, width=12, text="open sequence", command=self.clickopensequence)
        self.btn_opensequence.grid(column=0, row=0)
        
        self.btn_outputdir = Button(self.fileframe, width=12, text="output directory", command=self.clickoutputdir)
        self.btn_outputdir.grid(column=0, row=1)
        
        self.ent_opensequence = Entry(self.fileframe, width=40, textvariable=self.inputsequence_file)
        self.ent_opensequence.grid(column=1, row=0)
        
        self.ent_outputdir = Entry(self.fileframe, width=40, textvariable=self.outputdirectory_dir)
        self.ent_outputdir.grid(column=1, row=1)
        
        self.lbl_models = Label(self.fileframe, width=12, text="Models:")
        self.lbl_models.grid(column=0, row=2, sticky=W)
        
        self.cmb_models = ttk.Combobox(self.fileframe, state='readonly', textvariable=self.model)
        self.cmb_models.grid(column=1, row=2, sticky=W)
        self.cmb_models['values'] = self.populatemodels()
        
        self.btn_loadsettings = Button(self.fileframe, width=12,text="Load settings", command=self.loadsettings)
        self.btn_loadsettings.grid(column=3, row=0, padx=(20,10), sticky=E)
        
        self.btn_setsettings = Button(self.fileframe, width=12, text="Set Settings", command=self.setsettings)
        self.btn_setsettings.grid(column=3, row=1, padx=(20,10),sticky=E)
        
        
        #parameter frame
    def initparamframe(self, parent):
        self.paramframe = LabelFrame(parent, text="Parameters")
        self.paramframe.grid(column=0, row=4, sticky=NW)
        
        self.lbl_kmers = Label(self.paramframe, text="Kmers:")
        self.lbl_kmers.grid(column=0, row=0)
        
        self.ent_kmers = Entry(self.paramframe, width=5, textvariable=self.kmers)
        self.ent_kmers.grid(column=1, row=0)
        
        self.lbl_spacelen = Label(self.paramframe, text="Spacer length:")
        self.lbl_spacelen.grid(column=0, row=1)
        
        self.ent_spacelen = Entry(self.paramframe, width=5, textvariable=self.spacelen)
        self.ent_spacelen.grid(column=1, row=1)
        
        self.lbl_spacepos = Label(self.paramframe, text="Insert Spacer at:")
        self.lbl_spacepos.grid(column=0, row=2)
        
        self.ent_spacepos = Entry(self.paramframe, width=5, textvariable=self.spacepos)
        self.ent_spacepos.grid(column=1, row=2)
        
        self.lbl_startpos = Label(self.paramframe, text="Starting at")
        self.lbl_startpos.grid(column=0, row=3)
        
        self.ent_startpos = Entry(self.paramframe, width=5, textvariable=self.startpos)
        self.ent_startpos.grid(column=1, row=3)
         
        self.lbl_endpos = Label(self.paramframe, text="Ending at")
        self.lbl_endpos.grid(column=0, row=4)
        
        self.ent_endpos = Entry(self.paramframe, width=5, textvariable=self.endpos)
        self.ent_endpos.grid(column=1, row=4)
       
        #selection frame
    def initselectframe(self, parent):

        self.selectframe = LabelFrame(parent, text="Graphs")
        self.selectframe.grid(column=1, row=4, sticky=NW) 
        
        self.chk_creategraphs = Checkbutton(self.selectframe, text="Create Graphs", variable=self.creategraphs)
        self.chk_creategraphs.grid(column=0, row=0)
        
        self.lbl_forrvals = Label(self.selectframe, text="For R-Values:")
        self.lbl_forrvals.grid(column=0, row=1)
        
        self.lbl_toppercent = Label(self.selectframe, text="Top % ")
        self.lbl_toppercent.grid(column=0, row=2)
        
        self.ent_toppercent = Entry(self.selectframe, width=5, textvariable=self.topcut)
        self.ent_toppercent.grid(column=1, row=2)
        
        self.lbl_botpercent = Label(self.selectframe, text="Bottom % ")
        self.lbl_botpercent.grid(column=0, row=3)
        
        self.ent_botpercent = Entry(self.selectframe, width=5, textvariable=self.botcut)
        self.ent_botpercent.grid(column=1, row=3)
        
        #extra option frame
    def initoptionsframe(self, parent):
        self.optionsframe = LabelFrame(parent, text="Options")
        self.optionsframe.grid(column=2, row=4, sticky=NW)
        
        self.chk_keepseqs = Checkbutton(self.optionsframe, text="Keep all Sequences")
        self.chk_keepseqs.config(variable=self.keepseqs, onvalue=True, offvalue=False)
        self.chk_keepseqs.grid(column=0, row=0, sticky=W)
        
        self.chk_keeprvals = Checkbutton(self.optionsframe, text="Keep all R-Values")
        self.chk_keeprvals.config(variable=self.keeprvals, onvalue=True, offvalue=False)
        self.chk_keeprvals.grid(column=0, row=1, sticky=W)
        
        #message box frame
    def initmsgframe(self, parent):
        self.msgframe = Frame(parent)
        self.msgframe.grid(column=0, row=5, columnspan=4)
        
        self.msgbox = Text(self.msgframe, height=10, width=50)#, state=DISABLED)
        self.msgbox.pack(side=LEFT, fill=Y)
        
        self.scroll = Scrollbar(self.msgframe)
        self.scroll.pack(side=RIGHT, fill=Y)
        self.scroll.config(command=self.msgbox.yview)
        self.msgbox.config(yscrollcommand=self.scroll.set)
        
        #bottom progress start/stop bar
    def initbottombar(self, parent):
        self.bottombarframe = Frame(parent)
        self.bottombarframe.grid(column=0, row=6, columnspan=4)
        
        self.btn_start = Button(self.bottombarframe, width=10, text="Start", command=self.start)
        self.btn_start.pack(side=LEFT)
        
        self.lbl_progress = Label(self.bottombarframe, width=10, textvariable=self.progress)
        self.lbl_progress.pack(side=LEFT)
        
        self.btn_stop = Button(self.bottombarframe, width=10, text="Stop")
        self.btn_stop.pack(side=LEFT)
       
       
       
        #BUTTON METHODS#####################
    def start(self):
        self.msgbox.insert(END, "Starting...") 
        self.script_thread.start()
#        self.batch_process.start()
        root.update_idletasks()
        
    def stop(self):
        self.msgbox.insert(END, "Stop Process")
        self.script_thread.terminate()
#        self.batch_process.terminate()
        sys.exit()
        
    def run_script(self):
        process = subprocess.Popen([MainWindow.SCRIPT_FILE], shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        outmsg = process.stdout.read()

        self.msgbox.insert(END,process.returncode)

    def loadsettings(self):
        config = ConfigParser.ConfigParser()
        config.read(MainWindow.CONFIG_FILE)
        
        self.inputsequence_path.set(config.get('Files', 'input file'))
        self.inputsequence_file.set(os.path.basename(config.get('Files', 'input file')))
        self.outputdirectory_path.set(config.get('Files', 'output directory'))
        self.outputdirectory_dir.set(os.path.basename(config.get('Files', 'output directory')))
        self.cmb_models.set(os.path.basename(config.get('Files', 'model file')))
        
        self.kmers.set(config.getint('Params', 'kmers'))
        self.spacelen.set(config.getint('Params', 'spacer length'))
        self.spacepos.set(config.getint('Params', 'spacer after position'))
        self.startpos.set(config.getint('Params', 'starting at position'))
        self.endpos.set(config.getint('Params', 'ending at position'))
        
        self.creategraphs.set(config.getboolean('Graph', 'create graphs'))
        self.topcut.set(config.getint('Graph', 'top'))
        self.botcut.set(config.getint('Graph', 'bottom'))
        
        self.keepseqs.set(config.getboolean('Options', 'keep all sequences'))
        self.keeprvals.set(config.getboolean('Options', 'keep all rvalues'))
        
    def setsettings(self):
        config = ConfigParser.ConfigParser()
        config.read(MainWindow.CONFIG_FILE)
        
        config.set('Files', 'input file', self.inputsequence_path.get())
        config.set('Files', 'output directory', self.outputdirectory_path.get())
        config.set('Files', 'model file', MainWindow.MODEL_DIR + self.model.get())
        
        config.set('Params', 'kmers', self.kmers.get())
        config.set('Params', 'spacer length', self.spacelen.get())
        config.set('Params', 'spacer after position', self.spacepos.get())
        config.set('Params', 'starting at position', self.startpos.get())
        config.set('Params', 'ending at position', self.endpos.get())
        
        config.set('Graph', 'create graphs', self.creategraphs.get())
        config.set('Graph', 'top', self.topcut.get())
        config.set('Graph', 'bottom', self.botcut.get())
        
        config.set('Options', 'keep all sequences', self.keepseqs.get())
        config.set('Options', 'keep all rvalues', self.keeprvals.get())
        
        ks = os.path.basename(config.get('Options', 'sequences file'))
        rv = os.path.basename(config.get('Options', 'rvalues file'))
        gr = os.path.basename(config.get('Options', 'graph data file'))
        
        config.set('Options', 'sequences file', self.outputdirectory_dir.get() + '/' + ks)
        config.set('Options', 'rvalues file', self.outputdirectory_dir.get() + '/' + rv)
        config.set('Options', 'graph data file', self.outputdirectory_dir.get() + '/' + gr)
        
        with open(MainWindow.CONFIG_FILE, 'wb') as configfile:
            config.write(configfile)
        
    def clickopensequence(self):
        self.inputsequence_path.set(tkFileDialog.askopenfilename())
        self.inputsequence_file.set(os.path.basename(self.inputsequence_path.get()))
        
    def clickoutputdir(self):
        self.outputdirectory_path.set(tkFileDialog.askdirectory())
        self.outputdirectory_dir.set(os.path.basename(self.outputdirectory_path.get()))
   
   
        #OTHER METHODS
    def populatemodels(self):
        modellist=[]
        for files in os.listdir(MainWindow.MODEL_DIR):
            modellist.append(files)
        return modellist
    
    
#load up the interface
#this part can be part of another file importhing this one
root = Tk()
app = MainWindow(root)
root.mainloop()