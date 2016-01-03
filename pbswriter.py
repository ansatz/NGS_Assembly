# #!/usr/bin/env parallel
import types
from functools import wraps
from collections import Iterable
from collections import defaultdict,OrderedDict
from unidecorator import *
import subprocess, shlex
import os
from string import Template
from pexpect import pxssh
from extremeconfig import *
from pygraphviz import *
from copy import deepcopy

# -P procfile (change procfile)
# trinity_notebook/qsub.py

#dict model=(assemblytype=['trinity','star','bowtie'],
 #           scoretype=['N50','isoform'],
  #          genometype=['cow','deer'])
        
class Opt(object):
    #over all files
    GlobalRange = dict(kmer=[],trim=[],qualtreshold=[])
    Ranges = dict(kmer=[1,30],trim=[1,10],qualthreshold=[70,95])
    OptN = dict(kmer=int, trim=int, qualthreshold=int)

    def __init__(self,readpairs,numRounds=8,pbspipe=[]):
        self.readpairs      = readpairs #dict{readpair:mkdir_path}
        self.optround       = optround
        self.assemblytype   = assemblytype
        self.scoretype      = scoretype
        self.genometype     = genometype
        for pbs in pbspipe:
            self.Pipe = pbs.optpipe 
    def __repr__(self):
        return str(self.optround) 
    def html(self,minutes=10):
        ''' create graph dict(files=[],functions=[],score=[],time=[],cpu=[])
            create graph of sns lm model
            serve an html file of graph, update ever X minutes
        '''
        pass
        
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
        pass
    def randParam(self):
        pass
    def optround(self):
        pass



@universal
class PBS(object):
    def __init__(self, function_cmd, read_pair=None,\
                    QUEUENAME="cri", NODES=1, CPU=5, PBSFILE=None, WALL='24:00:00', JOBNAME=None, CALLDIR=None, OUTPUTDIR=None,\
                    LOGFILE=None, CMD=[],\
                    cowref='/export/home/gsingh6/cow/cow.gtf.gz', cowgtf=None,\
                    assemblytype=None,scoretype=None,genometype=None,\
                    output=[],optpipe=None):
        wraps(function_cmd)(self)
        #--pbs
        self.QUEUENAME   = QUEUENAME
        self.NODES       = NODES
        self.CPU         = CPU
        self.WALL        = WALL
        self.name=function_cmd.__name__
        if PBSFILE: 
            self.PBSFILE=PBSFILE
        else: 
            self.PBSFILE=self.name + 'pbs'
        if JOBNAME: 
            self.JOBNAME = JOBNAME
        else: 
            self.JOBNAME = self.name
        if LOGFILE:
            self.LOGFILE     = LOGFILE
        else:
            self.LOGFILE = self.name +'log'
        self.CALLDIR     = CALLDIR
        
        self.CMD            = CMD
        self.optpipe        = optpipe #pass this to Opt class, shallow copy
        self.Labels         = []
        self.graphdict      = None

        #--data
        self.OUTPUTDIR      = OUTPUTDIR
        self.output         = output
        self.DATA           = read_pair
        self.cow            = cowref         
        self.gtf            = cowgtf
        self.assemblytype   = assemblytype
        self.scoretype      = scoretype
        self.genometype     = genometype

        self.function_cmd = function_cmd
        self.params=dict(group=['denovo','align'], theta=['kmer','qual'])
        self.var= vars(self)

    def __call__(self,*args,**kwargs):
        args = OrderedDict(self.function_cmd(*args,**kwargs))
        self.graphdict = args
        print '***', self.graphdict
        for k,v in args.items():
            self.Labels.append(k)
            self.CMD.append(v)

    def once(self,cmd):
        '''remove cmd from Opt class, only present in PBS class
        '''
        self.optpipe= [ i for i in self.CMD if cmd == cmd]
        print 'once **'
        print len(self.optpipe)
        print len(self.CMD)
    
    def output(self,*args,**kwargs):
        '''append PBS.output at index n, mod string to output[n]
        '''
        if isinstance(args, Iterable):
            for s in args:
                self.output.append(s)
        elif n != None:
            self.output[n]=output


    def testremote(self):
        '''/TestCmds/126L_AGTCAA_L006_R1_001.fastq_head  126L_AGTCAA_L006_R1_001.fastq_head_1000 
        '''
        # - login
        s = pxssh.pxssh()
        user=config['user']; pw=config['pw']
        host1=config['host1']; host2=config['host2']
        host=host1

        # - writePBS (local)
        localpwd= os.getcwd()
        testtemp = localpwd+'/testtemp'
        self.writePBS(local=testtemp, test=True)
        pbsfile= testtemp+'/pbsfiles/' + str(self.PBSFILE)

        # - mkdir (server)
        testdir = '/export/home/gsingh6/TestCmds/' + self.name
        self.OUTPUTDIR = testdir
        cmd1 = 'mkdir ' + testdir + '/;' 
        cmd11= 'cd ' + testdir +';'

        # - shlex scp (local)
        scp = 'sshpass -p ' +pw+ ' scp '+pbsfile+' '+ user+'@'+host+':'+testdir 
        scplist = shlex.split(scp)
        #p = subprocess.Popen(scplist) 
        # - qsub (server)
        cmd2 = 'sleep 5;qsub '+ testdir+'/'+ str(self.PBSFILE) +';'
        cmd3 = 'qstat -a;'
         
        if not s.login(host, user, pw,auto_prompt_reset=False):
            print "SSH session failed on login."
            print str(s)
        else:
            print "SSH session login successful"
            s.prompt()
            s.sendline(cmd1)
            p = subprocess.Popen(scplist) 
            s.sendline(cmd11)
            s.sendline(cmd2 )
            s.sendline(cmd3 )
            s.prompt()         # match the prompt
            data = s.before     # print everything before the prompt.
            print data
            s.logout()

    def flatten(self,coll):
        for i in coll:
            if isinstance(i, Iterable) and not isinstance(i, basestring):
                for subc in self.flatten(i):
                    yield subc
            else:
                yield i 
    def writePBS(self, home='/export/home/gsingh6',local=None,test=False,opt=False):
        #if opt==True:
        #    commands=deepcopy(self.optpipe)
        #else:
        #    commands=deepcopy(self.CMD)

        _HASHBANG    = '#!/bin/bash'
        _JOBNAME     = '#PBS -N ' + str(self.JOBNAME)
        if test==True:
            _QUEUE       = '#PBS -q batch'
            #_OUTPUTDIR = '/export/home/gsingh6/TestCmds/' + self.name
            _OUTPUTDIR   = self.OUTPUTDIR
        else:
            _QUEUE       = '#PBS -q ' + str(self.QUEUENAME)
            _OUTPUTDIR   = self.OUTPUTDIR
            self.CMD.append('source /export/home/gsingh6/.bashrc')

        _NODES       = '#PBS -l nodes=' + str(self.NODES) + ':ppn=' + str(self.CPU)
        _WALL        = '#PBS -l walltime='+ str(self.WALL)
        _OUTPUT      = '#PBS -j oe' #log output and error
        if(self.CALLDIR): #WORKDIR
            self.CMD.append('cd ' + str(self.CALLDIR))
            _PBSFILE     = str(self.CALLDIR) + '/' + str(self.PBSFILE)
            _LOGFILE     = '#PBS -o ' + ''.join([self.CALLDIR, '/', self.LOGFILE ])
        else:
            if test==True:
                _PBSFILE = '/home/solver/NGS_Assembly/testtemp/pbsfiles/' + str(self.PBSFILE) #local
                _LOGFILE     = '#PBS -o ' + ''.join(['/export/home/gsingh6/TestCmds/',self.name,'/',self.LOGFILE ]) #server
            else:
                _PBSFILE     = home+'/pbsfiles/' + str(self.PBSFILE)
                _LOGFILE     = '#PBS -o ' + ''.join([home,'/pbslog/',self.LOGFILE ])
        print 'test**', self.CMD[:3]
        cflat=self.flatten(self.CMD)
        ctemp=[ Template(c).substitute(self.var) for c in cflat ]
        _CMD='\n'.join(ctemp)
        pbslist = [_HASHBANG,_JOBNAME,_QUEUE,_NODES,_WALL,_OUTPUT,_LOGFILE,_OUTPUTDIR,_CMD]
        if not _OUTPUTDIR:
            pbslist.pop(7)
        else:
            pbslist[7]=str(_OUTPUTDIR)
        pbs='\n'.join(pbslist)
        self.dotplot(self.name)
        print '** writePBS() called... '
        with open(_PBSFILE, 'w') as fp:
            fp.write(pbs)

    def dotplot(self,file):
        C=AGraph()
        addedge=lambda x,y: C.add_edge(x,y)
        lastnode=None
        flatlist=[]
        
        for k,v in self.graphdict.items():
            flatlist.append(k)
            if isinstance(v,basestring):
                flatlist.append(v)
            elif isinstance(v,list):
                for vv in v:
                    flatlist.append(vv)

        for i in xrange(0,len(flatlist)-1):
            C.add_edge( flatlist[i], flatlist[i+1] )

        for k,v in self.graphdict.items():
            nodek = C.get_node(k)
            nodek.attr['shape']='rect'
            #node.attr['pos']="%f,%f!"%()
            nodek.attr['fontcolor']='red'
        C.draw(file+'.png',prog='dot')


    def callPBS(self):
        commandQsub = 'qsub ' + str(self.CALLDIR) + '/' + str(self.PBSFILE)
        args = shlex.split(commandQsub)
        print(args[0:3])
        p = subprocess.Popen(args )

    def mkdir(self,optround):
        '''<genome><read_pair>_<assemblytype>_<scoretype>
        '''
        path = self.genometype + self.read_pair + '_' + self.assemblytpe + '_' + self.scoretype
        if not path:
            sys.mkdir(path)
        optpath = path+'/'+optround
        sys.mkdir(optpath)

    

# chainable decorator functions 
def P(function=None, originate=None,c=2,rerun=False):
    if function:
        return _P(function)
    else:
        def wrapper(function):
            return _P(function,c,originate,rerun)
        return wrapper
def _PBS(function=None, **kwargs):
    if function:
        return PBS(function)
    else:
        def wrapper(function):
            return _PBS(function,*kwargs)
        return wrapper
def Pipe(obj):
    def wrap(*args,**kwargs):
        print 'wrapped'
        obj()
        return obj(*args,**kwargs)
    return wrap
    
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
        for c in Pipe.share_fncs:
            print c.__name__
        return obj(*args,**kwargs)
    return wrap

#def once(*args,**kwargs):
#    def dec(obj):
#        wraps(obj)
#        def wrap(*args,**kwargs):
#            obj._once(*args,**kwargs)
#            return obj(*args,**kwargs)
#        return wrap
#    return dec

def output(PBS,*args,**kwargs):
    def wrap(*args,**kwargs):
        PBS._output(*args,n=None)
        return PBS(*args,**kwargs)
    return wrap



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



#class _Pipe(_PBS):
#    ''' 1 pipe for each unique group/param combination
#        groups[denovo|align] param[kmer|qual] => 4 
#        maintain monotonic parameter
#        utility: path(monotonic, intermediatefiles), log, jobs(parallel)
#    #decorator run functions
#        @pipe : rerun if update, mark state as complete, share intermediate files
#        @path : files, cmds to an output
#    '''
#    share_fncs=[]
#
#    def __init__(self,share_fnc=None,group=None,param=None,state=None,*args,**kwargs):
#        self.group=group
#        self.param= param
#        self.state=state
#        self.share=[]
#        self.share_fnc=share_fnc
#        PBS.__init__(self,*args,**kwargs)
#        #super(_Pipe,self).__init__()
# 
#    def __call__(self):
#        '''load all files
#            run for group (ie trinity or star)
#            dispatch Score
#        '''
#        print 'call Pipe'    
#        self.share_fncs.append(self.share_fnc)
#        print self.share_fncs
#        return self.share_fncs
#
#    def _makedir(self):
#        '''./<readpair>/<group>/<theta>/val
#        '''
#        groupparam=self.group+'_'+self.param
#        path='./outputdir/' + groupparam + '/' + PBS.params[readpair] + '/' + PBS.params[theta] +'_'+ PBS.params[val]
#        os.mkdir(path) 
#
#    def _pipe(self):
#        '''state, paths to output
#        '''
#
#    def _path(self):
#        pass
#
#    def _parallel(self):
#        '''run read-pairs in parallel as well as different groups,params
#        '''
