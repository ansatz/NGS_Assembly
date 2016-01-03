#!/usr/bin/parallel --shebang-wrap python
from pbswriter import PBS
from unidecorator import *
import itertools
import HTSeq

# downstream analysis steps

@universal
class Score(object):
    CMD=[]
    def __init__(self, function=None):
        #super(Score, self).__init__()
        self.coverageMax = None
        self.strucVarMax = None
        self.fastqreads = fastqreads
        self.range = None
        self.optN = None
        #self.OUTPUTDIR=PBS.OUTPUTDIR
        self.function=function

    def __call__(self,*args,**kwargs):
        CMD.append(self.function(*args,**kwargs))
        @classmethod
        def f():
            self.function(*args,**kwargs)

    def __repr__(self):
        '''print all functions'''
        cmd = ' '.join(CMD,)
        slf = ' '.join(['coverageN50', 'isoform_cufflinks'])
        print cmd + slf
            
       

    # - pipe functions
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

    def coverage_wig(input,type=None,output=None):
        ''' read coverage vector => wig graph
        '''
        if type==None:
            type=input[-3:].upper()
            print type
        intype={'SAM': HTSeq.SAM_Reader(input),
                'BAM': HTSeq.BAM_Reader(input)
                }
        align_file = intype[type]
        # SequenceWithQualities obj class; slots aln,iv
        # coverage
        # auto add chromosome vectors as needed
        # stranded two vectors per chrom
        coverage = HTSeq.GenomicArray( "auto", stranded=True, typecode='i' )
        for alngt in itertools.islice(align_file,10):
            if alngt.aligned:
                #print 'aligned'
                coverage[ alngt.iv ] += 1
        print 'COVER', coverage
        #wiggle graph .bed
        if output==None:
            oplus=type+'plus.wig'
            ominus=type+'minus.wig'
        else:
            oplus=output+'plus.wig'
            ominus=output+'minus.wig'
        coverage.write_bedgraph_file(oplus, '+')
        coverage.write_bedgraph_file( ominus, '-')


    def genes_reads(annot, align, type=None):
        ''' count read by gene
            exon overlapped by read
            GenomicArrayOfSets
            /export/home/gsingh6/cow/cow.gtf.gz'
        '''
        if type==None:
            type=align[-3:].upper()
            #print type
        intype={'SAM': HTSeq.SAM_Reader(align),
                'BAM': HTSeq.BAM_Reader(align)
                }
        align_file = intype[type]
        # test end included if mod3 not equal 0
        gtf_file=HTSeq.GFF_Reader(annot, end_included=True)
        # exons #objects of class GenomicFeature
        exons = HTSeq.GenomicArrayOfSets( "auto", stranded=True )
        for feature in itertools.islice(gtf_file,5000):
            if feature.type=='exon':
                exons[feature.iv]+=feature.name
        # ** interval search,range **
        interval = HTSeq.GenomicInterval( "chr1", 151377661, 151490311, "+" )
        #print list( exons[interval].steps() )
        # intersection set ; multiple gene/features/steps() to read 
        # counts dict init with gene names
        counts = {}
        for feature in gtf_file:
            if feature.type == "exon":
                counts[ feature.name ] = 0
        # list(iset)[0] from intersection set of exons(over interval) with SAM file
        for alnmt in itertools.islice(align_file,1000):
            if alnmt.aligned:
                iset = None
                for iv2, step_set in exons[ alnmt.iv ].steps():
                    if iset is None:
                        iset = step_set.copy()
                    else:
                        iset.intersection_update( step_set )
                if len( iset ) == 1:
                    counts[ list(iset)[0] ] += 1
        for name in sorted(counts.keys()):
            if counts[name]>0:
                print( name, counts[name])
        for k,v in counts.items():
            print k,v

    #htseq-count (counts per gene)
    #htseq-qa

    # - scoring functions
    def coverageN50(self,assemblyfile):  #,out=self.OUTPUTDIR+'/Trinity.fasta'):
        ''' median of cumsum contig length
        order the contigs by length
        cdf
        cdf(totalsum/2)
        create a annotation file from Trinity.fasta, as scaffold to coding region
        https://github.com/trinityrnaseq/trinityrnaseq/wiki/Transcriptome%20Contig%20Nx%20and%20ExN50%20stats
        '''
        lens = [(s.name,len(s)) for s in HTSeq.FastaReader(assemblyfile)]
        sl=sorted(lens, key=lambda x,y:y)
        index=0
        for i,s in enumerate(sl):
            sm=sum(sl[0:i])
            if sm> total:
                print '*',i,sum(sl[0:i])
                index=i
                break
        n50=sl[index]
        #can output a graph perhaps with binned sums on x-axis and y-axis is the names of contigs in that bin
        return n50

        
        
    def isoform_cufflinks(self):
        ''' isoform count cufflinks => fpkm (anova-test) <hits.sam>
            -u multireadcorrect
            -b bias correct (ref fasta)
            -p cpus 
        '''
        cufflinks = 'cufflinks -u -b ' + reffasta + ' -p ' + cpus + ' -o ' + ' -GTF ' + gtf
        
        
    def exons_DEXSeq(self):
        ''' count reads/exon =>
        '''
        pass
    def kmer_sailfish(self):
        ''' coverage without align in min using kmer
            => %kmer mapped
        '''
        pass
        




    # - pipeline optimization functions 
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
    
        



