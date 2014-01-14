#Written by Jia Li
#Florida State University, Department of Biology
#07-25-2012

#################################
#This class exists so that you don't have to
#call a bunch of configuration get()
#if a new configuration is added, just set the appropriate
#instance variable and add it on.
#################################




import ConfigParser

class configs:
    def __init__(self):
        cf = ConfigParser.ConfigParser()
        cf.read('config.cfg')
        
        self.svm_script = cf.get('Resources', 'svm script')
        self.gen_script = cf.get('Resources', 'gen script')
        self.batch_script = cf.get('Resources', 'batch script')
        self.select_script = cf.get('Resources', 'select script')
        self.creategff_script = cf.get('Resources', 'create gff script')
        
        self.input_file = cf.get('Files', 'input file')
        self.model = cf.get('Files', 'model file')
        
        self.kmers = cf.getint('Params', 'kmers')
        self.spacelen = cf.getint('Params', 'spacer length')
        self.spacepos = cf.getint('Params', 'spacer after position')
        self.startpos = cf.getint('Params', 'starting at position')
        self.endpos = cf.getint('Params', 'ending at position')
        
        self.creategraphs = cf.getboolean('Graph', 'create graphs')
        self.topcut = cf.getint('Graph', 'top')
        self.botcut = cf.getint('Graph', 'bottom')
        
        self.seqfile = cf.get('Options', 'sequences file') 
        self.rvalfile = cf.get('Options', 'rvalues file') 
        self.selectedfile= cf.get('Options', 'select data file')