import os
from glob import glob
import pandas as pd
from pyfasta import Fasta
from subprocess import Popen

def prepare(infofile, refname, output, minimapbin, samtoolsbin, picardbin, threadnumber):
    """
    :param infofile: a description file of details of each sample, example: sample_info.txt
    :param refname: a fasta format of each targeted sequence, exaple: Samples_gene.fa
    :param output: folder of temporary files
    :param minimapbin: minimap bin path
    :param samtoolsbin: samtools bin bath
    :param picardbin: picard bin path
    :return:
    """

    # Read in sample information file
    datainfo=pd.read_table(infofile,index_col="Index")
    outputname = os.path.join(output, 'minimap_run.sh')
    documentdir = os.path.dirname(os.path.abspath(infofile))
    indexedref = os.path.splitext(os.path.basename(refname))[0]+'.mmi'

    # Begin writing minimap2 command
    outio = open(outputname,"w")
    for idx in datainfo.index:
        genefile = documentdir + '/' + datainfo.loc[idx]['Gene_File'] + '.txt'
        genetable = pd.read_table(genefile, index_col="Index")
        fqname = documentdir+'/'+datainfo.loc[idx]['Sample']+'.fastq'
        bamfile = datainfo.loc[idx]['Sample'] + '.bam'
        g1start = genetable.loc[idx]['start_1'] - 10
        g3end = genetable.loc[idx]['end_3'] + 10
        gaplen = g3end - g1start
        print(minimapbin,' -t ', str(threadnumber), ' --rmq=yes -a ', indexedref, ' ', fqname, ' | ',
            picardbin,' SortSam -I /dev/stdin -O ', bamfile, ' -SO coordinate', sep='', file=outio)
        print(samtoolsbin,' index ',bamfile, file=outio)
    outio.close()
    print("Minimap2 command load!")

    ###run bwa mem
    minimapcmd="bash minimap_run.sh"
    runbwaalign = Popen(minimapcmd, shell=True, cwd=output)
    runbwaalign.communicate()
    return True