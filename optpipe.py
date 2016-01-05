from pbswriter import *
import glob
import re

class Opt(object):
    ''' Randomly selects file and theta[kmer, trim, qual_threshold] for nRounds.
        parameterRange{} is updated monotonically, selecting successive optN.
        Optimization is complete after nRounds are performed. 
        _Output_: 
            - heatmap[nRounds_x_3]
            - stats{score,runtime,cpu}
        _Input_:
            - group (genome, scorefnc)
            - files (readpairs, synthetic, refgenome)
    '''
    #class variables
    GlobalRange = dict(kmer=[],trim=[],qualtreshold=[])
    Ranges = dict(kmer=[1,30],trim=[1,10],qualthreshold=[70,95])
    OptN = dict(kmer=int, trim=int, qualthreshold=int)
    def __init__(self,reads=None,numRounds=4,score=None):
        self.reads          = reads #dict{readpair:mkdir_path}
        self.numRounds      = numRounds
        self.currRound      = 1
        self.range          = dict(kmer=[1,30],trim=[1,10],qualthreshold=[70,95])
        self.score          = score
        #self.assemblytype   = assemblytype
        #self.scoretype      = scoretype
        #self.genometype     = genometype

    #def __repr__(self):
    #    return str(self.optround) 

    def check(self):
        '''Checks hpcluster nodes queue.'''
        # get cpus 50% 

    def randFP(self, data='/mnt/store3/clustcriinterns/for_gobind'):
        '''Randomly select file, theta to run.
            20 readpairs deer
        '''
        #get readpairs, random select
          #{key pairs : values [files]}
          #126N_GTCCGC_L006_R1_001.fastq.gz #126N_GTCCGC_L006_R2_001.fastq.gz
        filenames=[filename for filename in os.listdir(data)]
        os.chdir(data)
        fastqfiles = [f for f in glob.glob('*.fastq.gz')]
        p = re.compile(r'^(?P<keyname>[a-zA-Z0-9]+\_[a-zA-Z0-9]+\_[a-zA-Z0-9])')
        filenames=defaultdict(list)
        for f in fastqfiles:
            if p.search(f):
                key = re.match(p,f).groups('keyname')
                filenames[key].append(f)
        #random files(l,r)
        fidx=random.randint(0,len(filenames)-1) 
        itr=0
        for k,v in filenames.items():        
            if itr==fidx:
                randreads=v
            itr+=1
        #random param
        pidx=random.randint(0,2)
        params=[kmer,trim,qualthresh]
        randparam=params[pidx]
        return (randreads,randparam) 
    
    def nR(self,minX,maxY):
        '''Returns optN based on quad(sum-series=max).'''
        #scipy.quadratic
        pass 
    
    def heatmap(self): pass
 
    def scoreoutput(self):
        fout[self.score]
    def html(self): pass

    HOME='/home/solver/NGS_Assembly/readpairs'
    def treeb(self,randfile=None,randparam=None,branch=None):
        '''Calls over randfiles, randparam for numRounds times.
            Recursion is not based on parameter opt range. But new instances will load files to continue branch recursion.
        '''
        #load from file
        if os.isfile('MONO_Scores.json'):
            self.load_monotonic()
        if os.isfile('OPTN_Parameters.json'):
            self.load_optn()
        #base_case
        f,p=self.randFP(data=HOME) #
        if branch is None:
            MONO_SCR[randparam].append(kc_score)
            return self.tree(randfile=f, randparam=p, branch=self.OPTN_BRANCH[p])
        #while(branch):
        while(self.currRound < self.numRounds):
            self.currRound+=1
            #run pbs
            pbsObj=PBS(randfile,randparam)
            pbsObj.wopr() #write, opt cpu, runtime, htmlstats
            self.numRounds+=1
            #eval score
            kc_score = self.get_detonate(pbsObject.kcoutput)
            f,p=self.randFP(data=HOME)
            #left
            if kc_score > branch[0]:
                branch=self.OPTN_BRANCH[randparam].pop()
                self.MONO_SCR[randparam].append(kc_score)
                return self.tree(randfile=f,randparam=p, branch=self.OPTN_BRANCH[p])
            #right
            if kc_score < branch[0]:
                MONO_SCR[randparam].append(kc_score) #sort at end
                self.OPTN_BRANCH[randparam]=range(branch[0],branch[1],1) #new range
                return self.tree(randfile=f, randparam=p, branch=self.OPTN_BRANCH[p])
       #write file 
        


#    def tree(self, branch=None,randFile=None,randParam=None,scoreoutfile=None):
#        '''Runs over random files and params(greedy) to create monotonic_params for runs. 
#           Every call takes randFile, randParam.
#            branch is tree interval either increment (R)n+optN or (L)n+1 based on score
#        '''
#        (randFiles,randparam)=self.randFP(DATADIR='/home/solver/NGS_Assembly/readpairs')
#        currScore=None
#        scoreoutfile='kc_'+n+'.txt'
#        #branch:
#        if branch is None:
#            branch = nR(self.range[randParam])
#        #base_case:
#        if branch:
#            if currScore is None:
#            currScore = ;
#        while(branch):
#            score = scored[theta(PBS.output)]
#            #left_increase:
#
#            if score > self.coverageMax:
#                self.coverageMax=score
#                left=branch[1:] #shifted index right
#                return tree(branch=left)
#            #right_decrease:
#            if score < self.coverageMax:
#                right=xrange( branch[0], branch[1] ) #increment of 1 between last and curr range
#                return tree(branch=right)
#  



    def optTreeG():
        '''yield from optTreeG
        '''
    def randParam(self):pass

    def optround(self): pass
