#!/usr/bin/parallel --shebang-wrap python
from pbswriter import PBS
from unidecorator import *

#@universal

#class Test(PBS):
#    def __init__(self,data='126L_AGTCAA_L006_R1_001.fastq_head_1000', pbs='126L', log='126LLog', home='.'):
#        self.data=data
#        self.pbs=pbs
#        self.log=log
#        self.home=home
#        #super(Test, self).__init__()
#        #self.test_fnc = test_fnc
#        #self.test_writePBS = 
#
#    def __call__(self,obj):
#        def wrap():
#            obj.DATA=self.data
#            obj.PBSFILE=self.pbs
#            obj.LOGFILE=self.log
#            obj.writePBS(self.home)
#        return wrap
#

@universal
class Score(object):
    def __init__(self,fastqreads,id):
        super(Score, self).__init__()
        self.coverageMax = None
        self.strucVarMax = None
        self.fastqreads = fastqreads
        self.range = None
        self.optN = None
    #run pip functions
    def qsubcmd(self, program, theta):
        '''return cmd string of PBS file
           parallel <program> <parameters>
        '''
        p = 'parallel --results' 
        #parallel --result self.fastq --kmer {KMER} ::: KMER {param} 

    def qsubpbs(self):
        ''' write PBS file
        '''
        pass
    #scoring functions
    #def coverageN50(self,out=PBS.OUTPUTDIR+'/Trinity.fasta'):
    #    pass
    def structVar(self):
        pass
    
    def optRange(self, range):
        ''' given range, returns optN size for worst-case optimization
            based on integer series sum n+(n+1)+...(n+i)=X, where optimal selection increases by 1 for each iteration of search
            opt range is (n^2+n)/2=X
            return list [1,5,10,15] for range(1,15)
        '''
        #range, worstcase is for [1,15] is 5
        #baseRange = [5,10,15] #TODO solve quadratic
        optN=5 #TODO solve quadratic
        optbranch = xrange(range[0],range[1],optN)
        return optbranch

    #recursive opt function
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
    #average over reads for a theta
    #possible thetas are kmer, trimmomatic, normalization

#make dir structure
    
        



