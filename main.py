from score import *
from pbswriter import *
from unidecorator import *
import os

@PBS
def twobit():
    dir_twoBit = '/export/home/gsingh6/cow/'
    twobit='twoBitToFa ' + dir_twoBit + 'bosTau8.2bit -bed=' + dir_twoBit + 'bosTau8.trf.bed.gz ' + dir_twoBit + 'bosTau8_exons.fa'
    p = deepcopy(paramsdict)
    p['pbsfilename']= p['job_name']= p['output_file'] = '2bittest'
    p['cpu']=5; p['CMD']=twobit; p['call_dir']=dir_twoBit
    writeTemplate(**p)

@PBS
def trimmomatic():
    ''' trimmomatic
    '''
    pass

@PBS
def bowtie():
    '''## bowtie (exons non-junction align)
     [burrows-wheeler index](https://www.youtube.com/watch?v=4WRANhDiSHM)
    '''
    a1='cd /export/home/gsingh6/$TRINITYOUTPUT'
    #samtools stats ../126L_AGTCAA_L006_R1_001.fastq > samtools_stats
    
    #cd /export/home/gsingh6/cow
    ##twoBitToFa bosTau8.2bit bosTau8_introns.fa
    
    ### indexdir='/export/home/gsingh6/cow/cow_index_introns'
    ### bowtie-build bosTau8_introns.fa "$indexDir" 
    
    ## bowtie2-build <fasta1, fasta2...> basename =>output basename.1.2.3.bt2 and rev.bt2
    ##cd /export/home/gsingh6/cow/cow_index_introns
    #cd /export/home/gsingh6/cow/cow2_idx_intron
    
    #echo -e "\e[31mbowtie2 building index <cow2_idx_intron>...\e[0m"
    #bowtie2-build ../bosTau8_introns.fa cow2_idx_intron
    #echo -e "\e[31mIndex <cow2_idx_intron> complete\e[0m"
    

@PBS
def tophat(th):
    pbit    ='/export/home/gsingh6/cow'
    tb      ='twoBitToFa '+ pbit + '/bosTau8.2bit bosTau8_introns.fa'
    idxdir  ='/export/home/gsingh6/tophat_test/cow_idx_intron'
    build   = 'bowtie2-build bosTau8_introns.fa' +  idxdir 
        
    datadir="/mnt/store3/clustcriinterns/for_gobind"
    R1="$datadir/126L_AGTCAA_L006_R1_001.fastq.gz"

    R2="$datadir/126L_AGTCAA_L006_R2_001.fastq.gz"
    ## run path/to/indexfile readfiles.fq
    th['CMD'].append('tophat' + idxdir + ' ' +  R1 + ' ' +  R2  + '--num-threads ' + th['CPU']) 

@PBS
def cufflinks():
    ''' isoform count cufflinks => fpkm (anova-test) <hits.sam>
        -u multireadcorrect
        -b bias correct (ref fasta)
        -p cpus 
    '''
    l = 'Cufflinks'
    c = 'cufflinks -u -b $cow -p $CPU -o $OUTPUTDIR -GTF $gtf'
    return {l:c}
cufflinks()
#cufflinks.OUTPUTDIR='cflnk'
cufflinks.writePBS(home=os.getcwd())

@PBS
def test():
    l='test'
    c="testc"
    return {l:c}
test()


@PBS
def detonate():
    '''http://deweylab.biostat.wisc.edu/detonate/vignette.html
        input: assembly
        output: prefix.score (
    '''
    L0='detonate'
    u0='cd /export/home/gsingh6/detonate-1.10'

    L1='prior refgenome mean,std len'
    u1='rsem-eval-estimate-transcript-length-distribution $cow $output'
    
    L2='rsem-eval'
    # --args [reads] [assembly] [output-prefix] [avg paired fragment len] [prior distribution params] [-p threads]
    c3='rsem-eval-calculate-score examples/toy_SE.fq examples/toy_assembly_1.fa examples/rsem_eval_1 76 --transcript-length-parameters rsem-eval/true_transcript_length_distribution/mouse.txt -p $CPU'
    
    L4="estimate 'true' assembly"
    # - rsem (not calc expr but calc alignment posterior prob)
    u5='rsem-prepare-reference --bowtie examples/toy_ref.fa examples/toy_rsem_ref'
    u6='rsem-calculate-expression -p $CPU examples/toy_SE.fq examples/toy_rsem_ref examples/toy_rsem_expr'
    # --args [--reference rsem output] [--expression rsem output] [--assembly output prefix] [--a-p-best | --a-p-sample]
    c7='ref-eval-estimate-true-assembly --reference examples/toy_rsem_ref --expression examples/toy_rsem_expr --assembly examples/ta --alignment-policy best'
    
    L8='kmer compression (kc) score'
    # - rsem estimate expr level of each sequence in estimated 'true' assembly
    u9='rsem-prepare-reference --bowtie examples/ta_0.fa examples/ta_0_ref'
    u10='rsem-calculate-expression -p 12 examples/toy_SE.fq examples/ta_0_ref examples/ta_0_expr'
    c11='ref-eval --scores kc --A-seqs examples/toy_assembly_1.fa --B-seqs examples/ta_0.fa --B-expr examples/ta_0_expr.isoforms.results --kmerlen 76 --readlen 76 --num-reads 46988 | tee examples/kc_1.txt'    
    
    L9='alignment based scores'
    # --contig and nucleotide F1 scores for align each assembly to the estimated "true" assembly, and vice versa, using Blat
    # [assembly] <=> [estimated-true] => output
    u13='blat -minIdentity=80 examples/ta_0.fa examples/toy_assembly_1.fa examples/toy_assembly_1_to_ta_0.psl'
    u14='blat -minIdentity=80 examples/toy_assembly_1.fa examples/ta_0.fa examples/ta_0_to_toy_assembly_1.psl'
    # -- now compute the contig and nucleotide scores
    c15='ref-eval --scores contig,nucl --weighted no --A-seqs examples/toy_assembly_1.fa --B-seqs examples/ta_0.fa --A-to-B examples/toy_assembly_1_to_ta_0.psl --B-to-A examples/ta_0_to_toy_assembly_1.psl --min-frac-identity 0.90 | tee examples/contig_nucl_1.txt' 
    
    #return {L0:u0, L1:[u1], L2:[c3], L4:[u5,u6],L8:[u9,u10,c11], L9:[u13,u14,c15]}
    return [(L0,u0), (L1,[u1]), (L2,[c3]), (L4,[u5,u6]), (L8,[u9,u10,c11]), (L9,[u13,u14,c15]) ]

detonate()  
#detonate.once('rsem-eval-estimate-transcript-length-distribution $cow $output[0]')

print '**main**'
#print detonate.optpipe
#print detonate.CMD[:3]
detonate.writePBS(home='/home/solver/NGS_Assembly')
#detonate.testremote()  

 
#http://informatics.fas.harvard.edu/rna-seq-data-analysis-2/
#
#@parallel
#@PBS
#def trinity():
#    '''
#    Trinity output
#    rsem
#    fpkm
#    
#    '''
#    putil='/export/share/apps/trinityrnaseq-2.0.2/util'
#    # global vars
#    source /export/home/gsingh6/.bashrc # modules load
#    R1='/export/home/gsingh6/126L_AGTCAA_L006_R1_001.fastq.gz'
#    R2='/export/home/gsingh6/126L_AGTCAA_L006_R2_001.fastq.gz'
#    putil='/export/share/apps/trinityrnaseq-2.0.2/util'
#    #output directories
#    trioutdir='/export/home/gsingh6/trinity_output'
#    rsemdir=$trioutdir'/RSEM_outdir'
#    blastdb='/export/home/gsingh6/blastdb'
#    trinotate='/export/home/gsingh6/trinotate_output'
#    transdecode='/export/home/gsingh6/TransDecoder/Trinity.fasta.transdecoder_dir'
#    # TRINITY
#    Trinity --seqType fq --left $R1 --right $R2 --CPU 30 --output $tridir --max_memory 200G
#    #sleep 2
#    ####rsem --dir RSEM_outdir bowtie_out
#    ###
#    ####../../util/align_and_estimate_abundance.pl --transcripts trinity_out_dir/Trinity.fasta --seqType fq --left reads.left.fq --right reads.right.fq --SS_lib_type RF --est_method RSEM --aln_method bowtie --trinity_mode --prep_reference --output_dir RSEM_outdir
#    ###
#    ###cd /export/home/gsingh6/trinity_output/
#    ###$putil/align_and_estimate_abundance.pl --transcripts /export/home/gsingh6/trinity_output/Trinity.fasta --seqType fq --left $R1 --right $R2 --est_method RSEM --aln_method bowtie2 --trinity_mode --prep_reference --output_dir $rsemdir
#    ###
#    ###echo '\n**rsem** complete\n'
#    ###
#    #### fpkm threshold plot
#    ###$putil/misc/count_features_given_MIN_FPKM_threshold.pl $rsemdir/RSEM.genes.results > cumul_counts.txt
#    ###
#    #### alignment stats  >bowtie_out
#    ###"$putil"/bowtie_PE_separate_then_join.pl  --seqType fq  --left $R1 --right $R2 --target $trioutdir'/Trinity.fasta' --aligner bowtie2
#    ###
#    ###"$putil"/SAM_nameSorted_to_uniq_count_stats.pl bowtie_out/bowtie_out.nameSorted.bam
#    ###
#    ####>2  bioconductor edgeR 
#    #### -- pool all run trinity single (unless bio replicates-> then indpt)
#    ####module load apps/R-3.2.0
#    ####./loadr.r
#    ###
#    # coding region; longest orf
#    # 6-mer (6open reading frames) HMM max_lkl(code|non-code) =>
#    #tdc='/export/share/apps/trinityrnaseq/trinity-plugins/transdecoder/TransDecoder'
#    
#    #../../trinity-plugins/transdecoder/TransDecoder -t Trinity.fasta -m 50
#    #cd /export/home/gsingh6/TransDecoder
#    #./TransDecoder.LongOrfs -t $trioutdir/Trinity.fasta
#    #TransDecoder -t Trinity.fasta -m 50
#    
#    echo See best_candidates.\*  for candidate ORFs


@PBS
def blasthomology():
    # blast homology search
    # search Trinity transcripts
    #blastx -query $trioutdir/Trinity.fasta -db $blastdb/uniprot_sprot.trinotate.pep -num_threads 20 -max_target_seqs 1 -outfmt 6 > $trinotate/blastx.outfmt6
    
    # search Transdecoder-transcripts
    #blastp -query $transdecode/longest_orfs.pep -db $blastdb/uniprot_sprot.trinotate.pep -num_threads 20 -max_target_seqs 1 -outfmt 6 > $trinotate/blastp.outfmt6
    
    ## uniref90 blast
    #blastx -query $trioutdir/Trinity.fasta -db $blastdb/uniprot_uniref90.trinotate.pep -num_threads 20 -max_target_seqs 1 -outfmt 6 > $trinotate/uniref90.blastx.outfmt6
    #blastp -query $transdecode/longest_orfs.pep -db $blastdb/uniprot_uniref90.trinotate.pep -num_threads 20 -max_target_seqs 1 -outfmt 6 > $trinotate/uniref90.blastp.outfmt6
    
    #cd $trinotate #uniprot.sprot.fasta analysis
    #  determine coverage length of tophits to target transcript
    #  args <blast+query> <makeblastdb -in> <blast+out>
    #  output files $trinotate/ blastx.outfmt6.hist  blastx.outfmt6.hist.list  blastx.outfmt6.w_pct_hit_length
    #blastx 
    #$putil/analyze_blastPlus_topHit_coverage.pl $trinotate/blastx.outfmt6 $trioutdir/Trinity.fasta $blastdb/uniprot_sprot.trinotate.pep
    #blastp
    #$putil/analyze_blastPlus_topHit_coverage.pl $trinotate/blastx.outfmt6 $trioutdir/Trinity.fasta $blastdb/uniprot_sprot.trinotate.pep
    
    # >characterize functional annotation
    #protein domain hmmer
    #hmmscan --cpu 20 --domtblout TrinotatePFAM.out $blastdb/Pfam-A.hmm $transdecode/longest_orfs.pep > pfam.log
    pass
    

#hmm search

@PBS
def debruijn():
    '''http://trinityrnaseq.sourceforge.net/advanced_trinity_guide.html#Butterfly_reconstruction
    '''
    pass

#@parallel
@PBS
def star():
    s='star'
    return [s]

@PBS
def khmer():
    '''meta transcriptome strip-split
        rna scaffolding rnapath
        http://ivory.idyll.org/blog/trinity-in-silico-normalize.html
     '''
    pass

#@Score
#def contigdistribution(lens,mapPercent):
#    pass
#
#@Score
#def homology(blastdb):
#    pass
#
#@Score
#def codingseq(hmmsec):
#    pass
#
#@Score(input='Trinity.fasta', output=['$sample.sam','$sample.bam','$sample_htseq_counts.txt'])

@PBS
def chado():
    '''http://angus.readthedocs.org/en/2014/drosophila_rnaseq_bwa_htseq.html
       create an annotation file that says that the entire length of each 'scaffold' is in fact a coding regio
    '''
    L1='chado'
    c1='cd $PBSDIR'
    c2='/mnt/ebs/tools/chado_test/chado/bin/gmod_fasta2gff3.pl \
            --fasta_dir $Trinityoutput \
            --gfffilename $Trinityoutput.gff3 \
            --type CDS \
            --nosequence'
    L2='map our paired-end sequence reads to the transcriptome reference'
    c3='bwa mem $reference /mnt/ebs/trimmed_x/$sample_1_pe /mnt/ebs/trimmed_x/$sample_2_pe > $sample.sam'
    c4='samtools view -Sb ${sample}.sam > ${sample}.unsorted.bam\
        samtools sort ${sample}.unsorted.bam ${sample}\
        samtools index ${sample}.bam'
    c5='htseq-count --format=bam --stranded=no --type=CDS --order=pos --idattr=Name ${sample}.bam Trinity_all_X.gff3 > ${sample}_htseq_counts.txt'
    return {L1:[c1], L2:[c2,c3,c4,c5]}


exit(0)
# --parallel over score, genome, assembly
Scoring = Score()
print Scoring.__repr__ # print score cmds, input{read_pairs Trinity.fasta} | output{gff sam etc}
Optpipe = Opt(pairs,Scoring)
print Optpipe.__repr__ # print rounds, graphviz_stats, 


Optpipe.run()
optpipe.test()
Optpipe.html(minutes=15)

'''
# --linear model
beta = Linear( Optpipe.stats )
beta.sns()
beta.stderr()
beta.anova()

# --debruijn
@PBS
def debruijn():
    # match()
    pass 

debruijn()
debruijn.graph()

originate='path to read_pairs'
optrinity=Opt(originate)
trinity()
print boost




# @write
#@Pipe
#@makedir
#@PBS
#def trinity():
#    a='trinity commands ...'
#    print a
#    return a
#trinity()
#print '**',trinity.CMD


#if __name__=="__main__":
#    ''' #paramsdict = dict(  pbsfilename=None, nodes= 1, cpu=5, job_name=None, que_name='cri', call_dir=None, output_file=None, CMD=None ); '''
#    # --runFastq('126L_AGTCAA_L006_R2_001.fastq.gz')
#
#    # --twoBit()
#    @call
#    @write
#    def twobit():
#        dir_twoBit = '/export/home/gsingh6/cow/'
#        twobit='twoBitToFa ' + dir_twoBit + 'bosTau8.2bit -bed=' + dir_twoBit + 'bosTau8.trf.bed.gz ' + dir_twoBit + 'bosTau8_exons.fa'
#        p = deepcopy(paramsdict)
#        p['pbsfilename']= p['job_name']= p['output_file'] = '2bittest'
#        p['cpu']=5; p['CMD']=twobit; p['call_dir']=dir_twoBit
#        writeTemplate(**p)
#
#    twobit()
#
#    # -- bowtie2
#    bt = deepcopy(paramsdict)
#    bt['call_dir'] = '/export/home/gsingh6/cow/cow_index_introns'
#    bt['CMD']= 'bowtie2-build ../bosTau8_introns.fa cow_index_introns'
#
#
#
#
#
#    # -- tophat
#
#
#    # -- trinity
#
     
def test_writePBS(obj):
    def wrap(*args,**kwargs):
        obj.DATA='126L_AGTCAA_L006_R1_001.fastq_head_1000'
        obj.PBSFILE='126Lf'
        obj.LOGFILE='126LLog'
        obj.writePBS(home='.')
        return obj(*args,**kwargs)
    return wrap

def test_callPBS(obj):
    def wrap(*args, **kwargs):
        print('test_call', obj)
        obj.callPBS()
        return obj(*args,**kwargs)
    return wrap 

#@test_callPBS
#@test_writePBS
#@Pipe
#@PBS
def testpbs():
    c='cmd 1'
    a='hola'
    return [c,a]
testpbs()
