from pbswriter import *

class Opt(object):
    ''' Randomly selects file and theta[kmer, trim, qual_threshold] for nRounds.
        parameterRange{} is updated monotonically, selecting successive optN.
        Optimization is complete after nRounds are performed. 
        _Output_: 
            - heatmap[nRounds,3]
            - stats{score,runtime,cpu}
        _Input_:
            - group (genome, scorefnc)
            - files (readpairs, synthetic, refgenome)
    '''
    #over all files
    GlobalRange = dict(kmer=[],trim=[],qualtreshold=[])
    Ranges = dict(kmer=[1,30],trim=[1,10],qualthreshold=[70,95])
    OptN = dict(kmer=int, trim=int, qualthreshold=int)

    def __init__(self,readpairs,numRounds=4):
        self.readpairs      = readpairs #dict{readpair:mkdir_path}
        self.numRounds       = numRounds
        self.optRound       = None
        self.assemblytype   = assemblytype
        self.scoretype      = scoretype
        self.genometype     = genometype
    def __repr__(self):
        return str(self.optround) 
        
    def optTree(self, branch=None, range=[1,15], theta='kmer'):
        '''return score value 
           theta is parameters, either kmer or quality
        '''
        scored = dict(kmer=self.coverageN50(),
                     quality = self.structVar())
        #branch:
        if branch==None:
            branch = ssOptRange(range)

        #base_case:
        if self.coverageMax==None:
            self.coverageMax = scored[theta(PBS.output)];
        
        while(branch):
            score = scored[theta(PBS.output)]
            #left_increase:
            if score > self.coverageMax:
                self.coverageMax=score
                left=branch[1:] #shifted index right
                return opttree(branch=left)
            #right_decrease:
            if score < self.coverageMax:
                right=xrange( branch[0], branch[1] ) #increment of 1 between last and curr range
                return opttree(branch=right)
   
    def optTreeG():
        '''yield from optTreeG
        '''
    def randParam(self):

    def optround(self):
