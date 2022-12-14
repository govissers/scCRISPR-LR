import os
from glob import glob
import pandas as pd
from pyfasta import Fasta
from subprocess import Popen

def prepare(infofile, refname, output, bwabin, samtoolsbin, picardbin, threadnumber):
    """
    :param infofile: a description file of details of each sample, example: sample_info.txt
    :param refname: a fasta format of each targeted sequence, exaple:Samples_gene.fa
    :param output: folder of temporary files
    :param bwabin: bwa bin path
    :param samtoolsbin: samtools bin bath
    :param picardbin: picard bin path
    :return:
    """

    # Read in sample information file
    datainfo=pd.read_table(infofile,index_col="Index")
    outputname = os.path.join(output, 'bwa_run.sh')
    documentdir = os.path.dirname(os.path.abspath(infofile))
    outio = open(outputname,"w")
    for idx in datainfo.index:
        genefile = documentdir + '/' + datainfo.loc[idx]['Gene_File'] + '.txt'
        genetable = pd.read_table(genefile, index_col="Index")
        fqname = documentdir+'/'+datainfo.loc[idx]['Sample']+'.fastq'
        bamfile = datainfo.loc[idx]['Sample'] + '.bam'
        g1start = genetable.loc[idx]['start_1'] - 10
        g3end = genetable.loc[idx]['end_3'] + 10
        gaplen = g3end - g1start
        print(bwabin,' bwasw -t ', str(threadnumber), ' -w ', gaplen, ' -q 0 ', os.path.basename(refname), ' ', fqname, ' | ',
            picardbin,' SortSam -I /dev/stdin -O ', bamfile, ' -SO coordinate', sep='', file=outio)
        print(samtoolsbin,' index ',bamfile, file=outio)
    outio.close()
    print("bwa command load!")

    ###run bwa mem
    bwacmd="bash bwa_run.sh"
    runbwaalign = Popen(bwacmd, shell=True, cwd=output)
    runbwaalign.communicate()
    return True