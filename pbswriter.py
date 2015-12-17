# #!/usr/bin/env parallel

from unidecorator import *
# -P procfile (change procfile)
# trinity_notebook/qsub.py
class MetaPBS(type):
    pass
@universal
class PBS():
    __metaclass__= MetaPBS
    def __init__(self, function_cmd, read_pair=None,\
                    QUEUENAME="cri", NODES=1, CPU=5, PBSFILE=None, WALL=None, JOBNAME=None, CALLDIR=None, OUTPUTDIR=None,\
                    LOGFILE=None, CMD=[]):
        self.QUEUENAME = QUEUENAME
        self.NODES = NODES
        self.CPU = CPU
        self.PBSFILE = PBSFILE
        self.WALL = WALL
        self.JOBNAME = JOBNAME
        self.CALLDIR = CALLDIR
        self.OUTPUTDIR = OUTPUTDIR
        self.LOGFILE = LOGFILE
        self.CMD = CMD
        self.DATA= read_pair
        
        self.function_cmd = function_cmd
        self.name=function_cmd.__name__
        self.params=dict(group=['denovo','align'], theta=['kmer','qual'])

    def __call__(self,*args,**kwargs):
        self.cmd = self.function_cmd(*args,**kwargs)
        #self.writeTemplate()

#    def writeTemplate2(self):
#        print 'writing your template...'
#        print 'self.name ', self.name
#        pass
#
#    def testP(self,readpair=None,theta=['kmer','quality'],group=['denovo','reference']):
#        ''' ./testoutputdir/<filepair>/<denovo|reference>/<kmer|quality>/
#            parallel for each readpair
#            for all possible
#        '''
#        theta=dict(kmer=self.monoKmer, quality=self.monoQuality)
#
#        c='parallel --results testoutputdir <command> --readpair {READPAIR} --kmer {KMER} --quality {QUALITY}\
#            ::: READPAIR ' + readpair + '::: KMER ' + kmer + '::: QUALITY ' + quality 
#        args = shlex.split(c)
#        p = subprocess.Popen(args )
       
    def writePBS(self, home='/export/home/gsingh6'):
        _HASHBANG    = '#!/bin/bash'
        _JOBNAME     = '#PBS -N ' + str(self.JOBNAME)
        _QUEUE       = '#PBS -q ' + str(self.QUEUENAME)
        _NODES       = '#PBS -l nodes=' + str(self.NODES) + ':ppn=' + str(self.CPU)
        _WALL        = '#PBS -l walltime='+ str(self.WALL)
        _OUTPUT      = '#PBS -j oe' #log output and error

        _CALLDIR     = self.CALLDIR
        _OUTPUTDIR   = self.OUTPUTDIR

        self.CMD.append('source ' + home + '/.bashrc')
        if(self.CALLDIR):
            self.CMD.append('cd ' + str(_CALLDIR))
            _PBSFILE     = str(_CALLDIR) + '/' + str(self.PBSFILE)
            _LOGFILE     = '#PBS -o ' + ''.join([_CALLDIR, '/', self.LOGFILE ])
        else:
            _PBSFILE     = home+'/pbsfiles/' + str(self.PBSFILE)
            _LOGFILE     = '#PBS -o ' + ''.join([home,'/pbslog/',self.LOGFILE ])

        _CMD = str(self.CMD)
        pbslist = [_HASHBANG,_JOBNAME,_QUEUE,_NODES,_WALL,_OUTPUT,_LOGFILE,_OUTPUTDIR,_CMD]
        if not _OUTPUTDIR:
            pbslist.pop(7)
        else:
            pbslist[7]=str(_OUTPUTDIR)
        pbs='\n'.join(pbslist)
        print '** writePBS'
        print self.DATA
        with open(_PBSFILE, 'w') as fp:
            fp.write(pbs)

    def callpbs(self):
        commandQsub = 'qsub ' + str(self.CALLDIR) + '/' + str(self.PBSFILE)
        args = shlex.split(commandQsub)
        p = subprocess.Popen(args )


 
#    def writeTemplate(self, pbsfilename, nodes, cpu, job_name, que_name, call_dir, output_file, CMD  ): 
#        ''' WRITES PBS FILE TEMPLATE WITH PARALLEL
#        '''
#        _PBSFILE     = '/export/home/gsingh6/pbsfiles/' + PBSFILE
#        _NODES       = 1
#        _HASHBANG    = '#!/bin/bash' 
#        _JOBNAME     = '#PBS -N ' + JOBNAME
#        _QUEUENAME   = '#PBS -q ' + QUENAME
#        _NODES       = '#PBS -l nodes=' + str(NODES) + ':ppn=' + str(CPU)
#        _WALL        = '#PBS -l walltime='+ str(WALL) 
#        _OUTPUT      = '#PBS -j oe' #log output and error
#        _LOGFILE     = '#PBS -o ' + ''.join(['export/home/gsingh6/pbslog/',LOGFILE]) 
#        _CALLDIR     = CALLDIR  
#        _OUTPUTDIR   = OUTPUTDIR
#        
#        if(call_dir): cd = 'cd ' + call_dir
#        modules = 'source /export/home/gsingh6/.bashrc' # \n cd $PBS_WORKDIR'
#        command = CMD
#       
#        #write pbs
#        pbs = '\n'.join([hashbang,jobname,queue,nodes,output,out,wall,cd,modules,command])
#        with open(pbsdir + pbsfilename, 'w') as fp:
#            fp.write(pbs)
#        #call pbs file 
#        #commandQsub = 'qsub -q ' + command
#        #commandQsub = 'qsub ' + pbsdir + pbsfilename
#        #print pbs
#        #args = shlex.split(commandQsub)
#        #p = subprocess.Popen(args )

class MetaPipe(type):
    pass
class MetaPP(MetaPBS,MetaPipe):
    pass
@universal
class Pipe():
    ''' 1 pipe for each unique group/param combination
        groups[denovo|align] param[kmer|qual] => 4 
        maintain monotonic parameter
        utility: path(monotonic, intermediatefiles), log, jobs(parallel)
    '''
    #__metaclass__=MetaPP
    share_fncs=[]

    def __init__(self,share_fnc=None,group=None,param=None,state=None):
        self.group=group
        self.param= param
        self.state=state
        self.share=[]
        self.share_fnc=share_fnc
    
 
    #decorator run functions
    ''' @compare : reruns w/o overwrite
        @rerun : overwrite
        @pipe : rerun if update, mark state as complete, share intermediate files
        @path : files, cmds to an output
    '''
    def __call__(self):
        self.share_fncs.append(self.share_fnc)
        return self.share_fncs

    def _makedir(self):
        '''./<readpair>/<group>/<theta>/val
        '''
        groupparam=self.group+'_'+self.param
        path='./outputdir/' + groupparam + '/' + PBS.params[readpair] + '/' + PBS.params[theta] +'_'+ PBS.params[val]
        os.mkdir(path) 

    def _pipe(self):
        '''state, paths to output
        '''

    def _path(self):
        pass

    def _parallel(self):
        '''run read-pairs in parallel as well as different groups,params
        '''
         

#decorator utility function(chainable to instance)
def makedir(obj):
    def wrap(*args,**kwargs):
        obj._makedir()
        return obj(*args,**kwargs)
    return wrap

def write(obj):
    #@wraps(obj)
    def wrap(*args,**kwargs):
        obj.writeTemplate2()
        print 'writing'
        return obj(*args,**kwargs)
    return wrap
  
def mark(obj):
    def wrap(*args,**kwargs):
        #obj.share_fncs.append(*args)
        for c in obj.share_fncs:
            print c.__name__
        return obj(*args,**kwargs)
    return wrap

            
# @write
#@Pipe
#@makedir
@PBS
def trinity():
    a='trinity commands ...'
    print a
    return a
trinity()
print '**',trinity

@mark
@Pipe
def boost():
    b='boost'
    return b

boost()
print trinity.name
