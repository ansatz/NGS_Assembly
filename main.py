from score import *
from pbswriter import *

#if __name__=="__main__":
#    '''paramsdict = dict(  pbsfilename=None, nodes= 1, cpu=5, job_name=None, que_name='cri', call_dir=None, output_file=None, CMD=None ); '''
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
        obj.callPBS()
        return obj(*args,**kwargs)
    return wrap 

@test_callPBS
@test_writePBS
@PBS
def testpbs():
    c='cmd 1'
    print c
    return c
testpbs()
