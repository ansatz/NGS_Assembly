#!/usr/bin/env parallel

from unidecorator import *
# -P procfile (change procfile)

@universal
class PBS(object):
    def __init__(self, dev_function_group, read_pair, theta,\
                cmd=1,originate=None,rerun=False,write=None):
        self.cmd = cmd
        self.dev_function_group = dev_function_group
        self.originate = originate
        self.name=dev_function_group.__name__
        self.rerun=rerun
        self.write=write
        self.output_coverage=None
        self.output_quality=None
        self.output_dir=None
        self.params=dict(readpair=None, group=['denovo','align'], theta=['kmer','qual'])

    def __call__(self,*args,**kwargs):
        self.cmd = self.dev_function(*args,**kwargs)
        #self.writeTemplate()

    def writeTemplate(self):
        print 'writing your template...'
        print 'self.name ', self.name
        pass

    def testP(self,readpair=None,theta=['kmer','quality'],group=['denovo','reference']):
        ''' ./testoutputdir/<filepair>/<denovo|reference>/<kmer|quality>/
            parallel for each readpair
            for all possible
        '''
        theta=dict(kmer=self.monoKmer, quality=self.monoQuality)

        c='parallel --results testoutputdir <command> --readpair {READPAIR} --kmer {KMER} --quality {QUALITY}\
            ::: READPAIR ' + readpair + '::: KMER ' + kmer + '::: QUALITY ' + quality 
        args = shlex.split(c)
        p = subprocess.Popen(args )
        
    def writeTemplate(self, pbsfilename, nodes, cpu, job_name, que_name, call_dir, output_file, CMD  ): 
        ''' WRITES PBS FILE TEMPLATE WITH PARALLEL
        '''
        _PBSFILE     = '/export/home/gsingh6/pbsfiles/' + PBSFILE
        _NODES       = 1
        _HASHBANG    = '#!/bin/bash' 
        _JOBNAME     = '#PBS -N ' + JOBNAME
        _QUEUENAME   = '#PBS -q ' + QUENAME
        _NODES       = '#PBS -l nodes=' + str(NODES) + ':ppn=' + str(CPU)
        _WALL        = '#PBS -l walltime='+ str(WALL) 
        _OUTPUT      = '#PBS -j oe' #log output and error
        _LOGFILE     = '#PBS -o ' + ''.join(['export/home/gsingh6/pbslog/',LOGFILE]) 
        _CALLDIR     = CALLDIR  
        _OUTPUTDIR   = OUTPUTDIR
        
        if(call_dir): cd = 'cd ' + call_dir
        modules = 'source /export/home/gsingh6/.bashrc' # \n cd $PBS_WORKDIR'
        command = CMD
       
        #write pbs
        pbs = '\n'.join([hashbang,jobname,queue,nodes,output,out,wall,cd,modules,command])
        with open(pbsdir + pbsfilename, 'w') as fp:
            fp.write(pbs)
        #call pbs file 
        #commandQsub = 'qsub -q ' + command
        #commandQsub = 'qsub ' + pbsdir + pbsfilename
        #print pbs
        #args = shlex.split(commandQsub)
        #p = subprocess.Popen(args )

@universal
class Pipe(PBS):
    ''' 1 pipe for each unique group/param combination
        groups[denovo|align] param[kmer|qual] => 4 
        maintain monotonic parameter
        utility: path(monotonic, intermediatefiles), log, jobs(parallel)
    '''

    def __init__(self,group=None,param=None,state=None):
        self.group=group
        self.param= param
        self.state=state
        self.share=[]
 
    #decorator run functions
    ''' @compare : reruns w/o overwrite
        @rerun : overwrite
        @pipe : rerun if update, mark state as complete, share intermediate files
        @path : files, cmds to an output
    '''
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
        obj.writeTemplate()
        print 'writing'
        return obj(*args,**kwargs)
    return wrap
  


            
@write
@PBS
def trinity():
    a='trinity commands ...'
    print a
    return a
trinity()
print trinity.cmd

