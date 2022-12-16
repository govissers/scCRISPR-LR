import os
from subprocess import Popen
from subprocess import PIPE
import re
import shutil
from CMlib import minimap2_run, subprocesspath
import time
import signal

# Test if bwa is working
def testminimap(minimapbin):
    """
    :param bwabin: bwa bin path
    :return: bool, True: bwa tested ok. False: bwa error
    """
    minimapcmd = [minimapbin]
    minimaprun = Popen(minimapcmd, stdout=PIPE, stderr=PIPE, shell=True)
    testres = False
    pat = re.compile('Version')

    for i in minimaprun.stderr.readlines():
        i = i.decode('utf-8').rstrip('\n')
        if re.search(pat, i):
            testres = True
    minimaprun.communicate()
    return testres


# Return the version of bwa used
def version(minimap):
    """
    :param minimap: bwa bin path
    :return: string, version of minimap
    """
    minimap = [minimap, '--version']
    minimap = Popen(minimap, stdout=PIPE, stderr=PIPE)
    i=minimap.stdout.readlines()[0]
    version = i.decode('utf-8').rstrip('\n')
    minimap.communicate()

    return version

# Index FASTA file of target genes
def index(minimapbin, reffile, samplefolder):
    """
    minimap index
    :param minimapbin: minimapbin bin path
    :param reffile: reference genome file
    :param samplefolder: sample dir
    :return: no retrun
    """
    refbasename = os.path.basename(os.path.abspath(reffile))
    refbasename_noext = os.path.splitext(refbasename)[0]
    dscopy = os.path.join(samplefolder, refbasename)
    shutil.copyfile(os.path.abspath(reffile), dscopy)
    
    minimapbin = os.path.abspath(minimapbin)
    minimapcmd = [minimapbin, '-d', refbasename_noext+'.mmi', refbasename]
    runminimapindex = Popen(minimapcmd, cwd=samplefolder)
    runminimapindex.communicate()


if __name__ == '__main__':


    bwapath = '../bin/bwa/x86_64-Darwin/bwa'

    seqlength = bwareflength(bwapath, '../Test/DM_404.fa')

    print(seqlength)

    # bwaalign(bwapath, '../Test/DM_404.fa', '../Test/Testsampe/oligo_tmp2.fa', '../Test/Testsampe/outfile.sam',4)
    # bwaindex(bwapath, '../Test/DM_404.fa', '../Test/Testsampe/')

    # bwapath = subprocesspath.subprocesspath(bwapath)
    #

    # seqlist = samfilter('../Test/Testsampe/outfile.sam', minas=45, maxxs=33)
    #
    # for i in seqlist:
    #
    #     print(i)

    # tester = bwaversion(bwapath)
    #
    # print(tester)
    #
    # res = bwaloci(bwapath, '../Test/Testsampe/DM_404.fa', '../Test/Testsampe/DM_test.faprobes.fa',threadnumber=4)
    #
    # for i in res:
    #     print(i)