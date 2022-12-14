import pysam
from pyfasta import Fasta
import matplotlib
from scipy import stats
from os import path
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import re
import os



def caldel(samfilename, start, end, genename):
    """

    :param samfilename: path and name of each bam file
    :param start: setting the start site of deletion calculation
    :param end: setting the end site of deletion calculation
    :param genename: the name of genome-editing target region
    :return:
    """
    n = 0

    mutateinfor = dict()

    deletelent = dict()
    samfile = pysam.AlignmentFile(samfilename, 'r')
    for read in samfile.fetch(genename):

        # print(read.cigartuples, read.cigarstring, read.reference_start, read.cigartuples[0][1], read.cigartuples[0][1]+read.reference_start)


        nowsite = read.reference_start
        #     print(read.cigarstring)
        for cigarnow in read.cigartuples:
            # print(cigarnow)
            cigartype = cigarnow[0]
            #         print(cigartype)
            cigarlenght = cigarnow[1]

            cigarend = nowsite + cigarlenght

            if start < nowsite < end:

                if cigartype == 2:

                    if cigarlenght < (end - start):

                        if cigarlenght in deletelent:

                            deletelent[cigarlenght] += 1
                        else:
                            deletelent[cigarlenght] = 1

            for i in range(nowsite, cigarend):

                if i in mutateinfor:

                    if cigartype in mutateinfor[i]:

                        mutateinfor[i][cigartype] += 1

                    else:

                        mutateinfor[i][cigartype] = 1

                else:

                    mutateinfor[i] = dict()

                    mutateinfor[i][cigartype] = 1

            nowsite += cigarlenght

        n += 1

    mutateinforpd = pd.DataFrame.from_dict(mutateinfor, orient='index').fillna(value=0)
    mutateinforpd['sum'] = mutateinforpd.sum(axis=1)
    mutateinforpd['delrate'] = mutateinforpd[2]/mutateinforpd['sum']
    deletelentpd = pd.DataFrame.from_dict(deletelent, orient='index')

    return (mutateinforpd, deletelentpd)

def plotpdf(groupinfo, refname, output, bamdir):

    """

    :param groupinfo: a description file of details of each group, example: group_infor.txt
    :param refname: a fasta format of the sequence in the target region, exaple:Samples_gene.fa
    :param output: folder of final result
    :param bamdir: folder of temporary files
    :return:
    """

    groupinfor = pd.read_table(groupinfo)
    groupinfor = groupinfor.dropna(axis=0, how='any')
    fa = Fasta(refname)

    for idx in groupinfor.index:

        repbam1 = os.path.join(bamdir, groupinfor.loc[idx]['rep1'] + '.bam')
        repbam2 = os.path.join(bamdir, groupinfor.loc[idx]['rep2'] + '.bam')
        ckbam = os.path.join(bamdir, groupinfor.loc[idx]['control'] + '.bam')
        strand = groupinfor.loc[idx]['strand']
        type = ''
        if (re.search("gRNA", groupinfor.loc[idx]['group'])):
            if strand == '+':
                start = groupinfor.loc[idx]['start'] - 10
                end = groupinfor.loc[idx]['end'] + 10
                type = 'gf'
            else:
                start = groupinfor.loc[idx]['start'] - 10
                end = groupinfor.loc[idx]['end'] + 10
                type = 'gr'

        elif (re.search("crRNA", groupinfor.loc[idx]['group'])):
            if strand == '+':
                start = groupinfor.loc[idx]['start']
                end = groupinfor.loc[idx]['end'] + 30
                type = 'cf'
            else:
                start = groupinfor.loc[idx]['start'] - 30
                end = groupinfor.loc[idx]['end']
                type = 'cr'
        genename = groupinfor.loc[idx]['gene']
        namenow = groupinfor.loc[idx]['group']


        if (path.exists(repbam1) and path.exists(repbam2)) and path.exists(ckbam):

            #print(repbam1, repbam2, ckbam, start, end, genename, namenow)

            seq = fa[genename][start - 1:end].upper()

            seqlist = list()

            for nt in seq:
                seqlist.append(nt)

            (mutateinforpd1, deletelentpd1) = caldel(samfilename=repbam1, start=start, end=end, genename=genename)
            (mutateinforpd2, deletelentpd2) = caldel(samfilename=repbam2, start=start, end=end, genename=genename)

            rep1 = mutateinforpd1.loc[start:end].delrate
            rep2 = mutateinforpd2.loc[start:end].delrate
            reg = pd.concat([rep1, rep2], axis=1)

            regmean = reg.mean(axis=1)
            stdrr = reg.sem(axis=1)

            (mutateinforpdCK, deletelentpdCK) = caldel(samfilename=ckbam, start=start, end=end, genename=genename)

            seqlistPAM = list()
            seqlistother = list()
            regPAM = list()
            if type == 'gf':
                seqlistPAM = seqlist[-13:-10]
                seqlistother = seqlist
                seqlistother[-13:-10] = ['', '', '']
                regPAM = regmean[-13:-10]
            if type == 'gr':
                seqlistPAM = seqlist[10:13]
                seqlistother = seqlist
                seqlistother[10:13] = ['', '', '']
                regPAM = regmean[10:13]
            if type == 'cf':
                seqlistPAM = seqlist[0:4]
                seqlistother = seqlist
                seqlistother[0:4] = ['', '', '', '']
                regPAM = regmean[0:4]
            if type == 'cr':
                seqlistPAM = seqlist[-4:]
                seqlistother = seqlist
                seqlistother[-4:] = ['', '', '', '']
                regPAM = regmean[-4:]


            fig, (ax0, ax1) = plt.subplots(ncols=2, sharey=True, figsize=(16, 9))
            # ax0.bar(regmean.index, regmean, yerr=stdrr)
            ax0.bar(regmean.index, regmean, color='purple')
            # add errorbar, elinewidth：errorbar line with; capsize/capthick:上下横线长短／粗细，ls:linestyle='None'去掉连接线。 ecolor: errorbar line color
            ax0.errorbar(regmean.index, regmean, yerr=stdrr, fmt='', elinewidth=0.5, capsize=2, capthick=0.5, ls='None',
                         ecolor='black')
            ax0.set_title(namenow)

            ax0.set_xticks(regmean.index, minor=True)
            ax0.set_xticklabels(seqlistother, color="black", minor=True)  # minor=True表示次坐标轴
            ax0.set_xticks(regPAM.index)
            ax0.set_xticklabels(seqlistPAM, color="red")

            # ax0.set_xticks(regmean.index)
            # ax0.set_xticklabels(seqlist)
            # ax0.tick_params(labelsize=8)


            ckname = namenow + ' Control'
            pdfname = os.path.join(output, namenow + '.pdf')
            regck = mutateinforpdCK.loc[start:end]
            ax1.bar(regck.index, regck.delrate, color='grey')
            ax1.set_title(ckname)
            ax1.set_xticks(regck.index,minor=True)
            ax1.set_xticklabels(seqlist,minor=True)
            ax1.set_xticks(regPAM.index)
            ax1.set_xticklabels(seqlistPAM, color="red")
            #ax1.tick_params(labelsize=8)
            #     plt.show()

            plt.savefig(pdfname)
            plt.close(fig)
            print("group",namenow, "finished!")
