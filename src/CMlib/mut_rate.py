from distutils.log import error
import os
from pickle import TRUE
import pysam
from pyfasta import Fasta
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import re
from glob import glob

def rate_cal(infofile, refname, output, bamdir):
    """
    :param infofile: a description file of details of each sample, example: sample_info.txt
    :param refname: a fasta format of the sequence in the target region, exaple:Samples_gene.fa
    :param output: folder of final result
    :param bamdir: folder of temporary files
    :return:
    """

    # Read in sample information
    info=pd.read_table(infofile,index_col="Index")

    # groupinfor = pd.read_table(groupinfo)
    # groupinfor.loc[:,pd.isnull(groupinfor).all()] = "UNKNOWN"
    # groupinfor=groupinfor.fillna("UNKNOWN")

    # Use group info to determine strandedness of read
    # stranddict = dict()
    # for idy in groupinfor.index:
    #     stranddict[groupinfor.loc[idy].rep1] = groupinfor.loc[idy].strand
    #     stranddict[groupinfor.loc[idy].rep2] = groupinfor.loc[idy].strand
    #     stranddict[groupinfor.loc[idy].control] = groupinfor.loc[idy].strand

    fa = Fasta(refname)

    print("Start mutation calculation!")

    for idx in info.index:

        # Get path of bam file
        bamname =  os.path.join(bamdir, info.loc[idx].Sample+'.bam')
        print("Calculating",bamname)

        # Get path of gene file
        genefilename = info.loc[idx].Gene_File
        infodir = os.path.dirname(infofile)
        genefile = os.path.join(infodir, genefilename+".txt")

        geneinfo = pd.read_table(genefile, index_col="Index")
        genelist = sorted(geneinfo.geneid)
        geneinfo = geneinfo.set_index("geneid")

        cas = info.loc[idx].Cas9
        sample = info.loc[idx].Sample

        # Read in BAM alignment file
        samfile = pysam.AlignmentFile(bamname, "rb", check_sq=False)

        # Generate list of all unique cells found in BAM alignment file
        cellids = []
        for read in samfile.fetch():
            cellid = read.query_name[0:8]
            if cellid not in cellids:
                cellids.append(cellid)

        ## Initialize cell-by-gene dfs
        # Mutation rate
        mtdf = pd.DataFrame(columns=cellids)
        mtdf.insert(0, 'genename', value=None)

        # Frameshift rate
        fsdf = pd.DataFrame(columns=cellids)
        fsdf.insert(0, 'genename', value=None)

        # Replacement rate
        repdf = pd.DataFrame(columns=cellids)
        repdf.insert(0, 'genename', value=None)

        # Insertion rate
        insdf = pd.DataFrame(columns=cellids)
        insdf.insert(0, 'genename', value=None)

        # Deletion rate
        deldf = pd.DataFrame(columns=cellids)
        deldf.insert(0, 'genename', value=None)

        # Insertion and Deletion rate
        insdeldf = pd.DataFrame(columns=cellids)
        insdeldf.insert(0, 'genename', value=None)

        # Large deletion rate
        largedeldf = pd.DataFrame(columns=cellids)
        largedeldf.insert(0, 'genename', value=None)

        # Loop through each gene in found in the FASTA file to generate rates of
        # each type of mutation for each cell
        for gene in genelist:

            print("Gene: " + gene)
            
            # Determine Cas9 editing range for each guide
            (start1, end1) = get_range(info, idx, geneinfo, gene, "1")
            (start2, end2) = get_range(info, idx, geneinfo, gene, "2")
            (start3, end3) = get_range(info, idx, geneinfo, gene, "3")
            
            cells = pd.DataFrame(
                {'cellids': cellids,
                'totalcov': [set() for _ in range(len(cellids))],
                'mtreads': [set() for _ in range(len(cellids))],
                'replacement': [set() for _ in range(len(cellids))],
                'insertion': [set() for _ in range(len(cellids))],
                'deletion': [set() for _ in range(len(cellids))],
                'insertion_deletion': [set() for _ in range(len(cellids))],
                'insertion_only': [set() for _ in range(len(cellids))],
                'deletion_only': [set() for _ in range(len(cellids))],
                'len_change': [dict() for _ in range(len(cellids))],
                'large_del': [set() for _ in range(len(cellids))]}
            )

            cells = cells.set_index('cellids')

            pileup = samfile.pileup(gene, max_depth=50000)
            for pileupcolumn in pileup:
                for pileupread in pileupcolumn.pileups:

                    if (start1 <= pileupread.query_position_or_next <= end1
                        or start2 <= pileupread.query_position_or_next <= end2
                        or start3 <= pileupread.query_position_or_next <= end3):

                        cellid = pileupread.alignment.query_name[0:8]
                        cells.loc[cellid]['totalcov'].add(pileupread.alignment.query_name)

                        # If neither deletion nor insertion, check reference for replacement
                        if not pileupread.is_del and not pileupread.is_refskip:
                            querybase = pileupread.alignment.query_sequence[pileupread.query_position]
                            refbase = fa[gene][pileupcolumn.pos].upper()

                            if querybase != refbase:
                                #                         replace += 1
                                cells.loc[cellid]['mtreads'].add(pileupread.alignment.query_name)
                                cells.loc[cellid]['replacement'].add(pileupread.alignment.query_name)

                        # Indel > 0 indicates insertion
                        if pileupread.indel > 0:
                            #                     insert += 1
                            cells.loc[cellid]['mtreads'].add(pileupread.alignment.query_name)
                            cells.loc[cellid]['insertion'].add(pileupread.alignment.query_name)
                            cells.loc[cellid]['len_change'][pileupread.alignment.query_name] = cells.loc[
                                cellid]['len_change'].get(pileupread.alignment.query_name, 0) + 1

                        # Indel < 0 indicates deletion
                        if pileupread.indel < 0:
                            #                     deletion += 1
                            cells.loc[cellid]['deletion'].add(pileupread.alignment.query_name)
                            cells.loc[cellid]['mtreads'].add(pileupread.alignment.query_name)
                            cells.loc[cellid]['len_change'][pileupread.alignment.query_name] = cells.loc[
                                cellid]['len_change'].get(pileupread.alignment.query_name, 0) - 1

            # Calclate the percentage of large deletions
            dels = [gene+'_DEL1', gene+'_DEL2', gene+'_DEL3']
            for deletion in dels:
                pileup = samfile.pileup(deletion, max_depth=50000)
                for pileupcolumn in pileup:
                    for pileupread in pileupcolumn.pileups:
                        cellid = pileupread.alignment.query_name[0:8]
                        cells.loc[cellid]['totalcov'].add(pileupread.alignment.query_name)
                        cells.loc[cellid]['large_del'].add(pileupread.alignment.query_name)

            # Calculate rates of each type of mutation
            cells = calc_mut_rates(cells)

            # Populate dataframes for each mutation type
            mtdf.loc[len(mtdf.index)] = fillrow(cells, 'percent_mt', gene)
            fsdf.loc[len(fsdf.index)] = fillrow(cells, 'percent_frameshift', gene)
            repdf.loc[len(repdf.index)] = fillrow(cells, 'percent_replacement', gene)
            insdf.loc[len(insdf.index)] = fillrow(cells, 'percent_insertion', gene)
            deldf.loc[len(deldf.index)] = fillrow(cells, 'percent_deletion', gene)
            insdeldf.loc[len(insdeldf.index)] = fillrow(cells, 'percent_insertion_deletion', gene)
            largedeldf.loc[len(largedeldf.index)] = fillrow(cells, 'percent_largedel', gene)


        # Make folder in output directory for given sample if it doesn't already exist
        if not os.path.isdir(os.path.join(output, sample)):
            os.mkdir(os.path.join(output, sample))

        # Save dataframe of each mutation type to results folder
        mtdfpath = os.path.join(output, sample, 'mt_rate.csv')
        mtdf.to_csv(path_or_buf=mtdfpath)

        fsdfpath = os.path.join(output, sample, "frameshift_rate.csv")
        fsdf.to_csv(path_or_buf=fsdfpath)

        repdfpath = os.path.join(output, sample, "replacement_rate.csv")
        repdf.to_csv(path_or_buf=repdfpath)

        insdfpath = os.path.join(output, sample, "insertion_rate.csv")
        insdf.to_csv(path_or_buf=insdfpath)
        
        deldfpath = os.path.join(output, sample, "deletion_rate.csv")
        deldf.to_csv(path_or_buf=deldfpath)

        insdeldfpath = os.path.join(output, sample, "insertion_deletion_rate.csv")
        insdeldf.to_csv(path_or_buf=insdeldfpath)

        lgdeldfpath = os.path.join(output, sample, "large_deletion_rate.csv")
        largedeldf.to_csv(path_or_buf=lgdeldfpath)
        
    samfile.close()

# Helper function to calcualate reads with frameshifts due to mutations in the specified region
# surrounding the PAM
def frameshift(len_change):
    """
    :param len_change: dictionary of reads and their respective differences in length from the reference sequence
            due to mutations
    :return: set() of reads with frameshifts
    """
    frameshifts = set()
    for read in len_change.keys():
        if not len_change[read] % 3 == 0:
            frameshifts.add(read)

    return frameshifts

# Helper function that populates a row from the cell dataframe into a new df for each mutation type
def fillrow(celldf, colname, genename):
    """
    :param celldf: DataFrame that contains the names of reads for each mutation type, for each cell
    :param colname: String name of the column from celldf to populate the row with values
    :param genename: String name of the gene corresponding to that celldf
    :return: list() containing 
    """
    row = celldf[colname].transpose().values
    row = row.tolist()
    row.insert(0, genename)

    return row


# Helper function that takes the dataframe of read sets to calculate rates of each
# mutation type
def calc_mut_rates(cells):
    """
    :param cells: DataFrame that contains the names of reads for each mutation type, for each cell
    :return: cells DataFrame with columns new containing rates of each mutation type
    """
    cells['insertion_deletion'] = cells.apply(lambda x: x.insertion.union(x.deletion), axis=1)
    cells['deletion_only'] = cells.apply(lambda x: x.deletion.difference(x.insertion), axis=1)
    cells['insertion_only'] = cells.apply(lambda x: x.insertion.difference(x.deletion), axis=1)
    cells['frameshift'] = cells['len_change'].apply(frameshift, 1)
    cells['percent_mt'] = cells.apply(lambda x: (len(x.mtreads) / len(x.totalcov)) , axis=1)
    cells['percent_replacement'] = cells.apply(lambda x: (len(x.replacement) / len(x.totalcov)) , axis=1)
    cells['percent_insertion'] = cells.apply(lambda x: (len(x.insertion) / len(x.totalcov)) , axis=1)
    cells['percent_deletion'] = cells.apply(lambda x: (len(x.deletion) / len(x.totalcov)) , axis=1)
    cells['percent_insertion_deletion'] = cells.apply(lambda x: (len(x.insertion_deletion) / len(x.totalcov)) , axis=1)
    cells['percent_deletion_only'] = cells.apply(lambda x: (len(x.deletion_only) / len(x.totalcov)) , axis=1)
    cells['percent_insertion_only'] = cells.apply(lambda x: (len(x.insertion_only) / len(x.totalcov)) , axis=1)
    cells['percent_frameshift'] = cells.apply(lambda x: (len(x.frameshift) / len(x.totalcov)) , axis=1)
    cells['percent_largedel'] = cells.apply(lambda x: (len(x.large_del) / len(x.totalcov)), axis = 1)
    
    return cells

# Takes a sample info dataframe and its assocaited gene info dataframe, as well as
# a gene and guide number of interest to calculate the range of possible CRISPR
# edits
def get_range(info, idx, geneinfo, gene, guidenum):
    """
    :param info: Sample information df
    :param idx: Current index of information df
    :param geneinfo: Gene information df
    :param gene: Current gene targeted by guides
    :param guidenum: Guide number to find range for
    :return: range of basepairs to search for edits
    """
    # Determine strandedness of guide
    strand = geneinfo.loc[gene]['direction_'+guidenum]

    if (re.search("gRNA", info.loc[idx].Cas9)):
            start = geneinfo.loc[gene]['start_'+guidenum] - 10
            end = geneinfo.loc[gene]['end_'+guidenum] + 10

    elif (re.search("crRNA", info.loc[gene].Cas9)):
        if strand == '+':
            start = geneinfo.loc[gene]['start_'+guidenum]
            end = geneinfo.loc[gene]['end_'+guidenum] + 30
        else:
            start = geneinfo.loc[gene]['start_'+guidenum] - 30
            end = geneinfo.loc[gene]['end_'+guidenum]

    return (start, end)
            

def display(groupinfo, output):
    """
    :param groupinfo: a description file of details of each group, example: group_infor.txt
    :param output: folder of final result
    :return:
    """

    mutfile = os.path.join(output, 'mut_rate.txt')
    mut_rate = pd.read_table(mutfile, sep='\t')
    groupinfor = pd.read_table(groupinfo)
    groupinfor = groupinfor.dropna(axis=0, how='any',thresh=6)  ##过滤表哥中没填满的行，thresh=7表示至少7个数不是NA,控制treatment和CK至少有一个
    #groupinfor.ix[:,pd.isnull(groupinfor).all()] = "UNKNOWN"
    groupinfor=groupinfor.fillna("UNKNOWN") ##填充表格中NaN处
    #print(groupinfor)
    mut_result = dict()
    for idx in mut_rate.index:
        mut_result[mut_rate.loc[idx].Sample] = mut_rate.values[idx]  ##读入mutation信息
    #    mut_result['OsPDS-RZ-gRNA1_Rep1'][2]

    ## prepare for display
    #replace = list()
    #replace_yerr = list()
    mutation=list()
    mutation_yerr = list()
    insertO = list()
    insertO_yerr = list()
    deletionO = list()
    deletionO_yerr = list()
    insert_deletion = list()
    insert_deletion_yerr = list()
    glist = list()
    ck_glist = list()

    ck_mutation = list()
    ck_insertO = list()
    ck_deletionO = list()
    ck_insert_deletion = list()
    for idy in groupinfor.index:
        rep1 = groupinfor.loc[idy].rep1
        rep2 = groupinfor.loc[idy].rep2
        ck = groupinfor.loc[idy].control
        if (mut_result.__contains__(rep1) and mut_result.__contains__(rep2)):

            # replace_mean = np.mean([mut_result[rep1][2], mut_result[rep2][2]])  ##np.mean([1,2,3,4,5])
            # #    print(group_mean)
            # replace.append(replace_mean)
            # replace_std = np.std([mut_result[rep1][2], mut_result[rep2][2]])  ## 标准差
            # #    print("std", group_var)
            # replace_yerr.append(replace_std)
            mutation_mean = np.mean([mut_result[rep1][1], mut_result[rep2][1]])  ##np.mean([1,2,3,4,5])
            #    print(group_mean)
            mutation.append(mutation_mean)
            mutation_std = np.std([mut_result[rep1][1], mut_result[rep2][1]])  ## 标准差
            #    print("std", group_var)
            mutation_yerr.append(mutation_std)

            insertO_mean = np.mean([mut_result[rep1][3], mut_result[rep2][3]])
            insertO.append(insertO_mean)
            insertO_std = np.std([mut_result[rep1][3], mut_result[rep2][3]])
            insertO_yerr.append(insertO_std)

            deletionO_mean = np.mean([mut_result[rep1][4], mut_result[rep2][4]])
            deletionO.append(deletionO_mean)
            deletionO_std = np.std([mut_result[rep1][4], mut_result[rep2][4]])
            deletionO_yerr.append(deletionO_std)

            insert_deletion_mean = np.mean([mut_result[rep1][5], mut_result[rep2][5]])
            insert_deletion.append(insert_deletion_mean)
            insert_deletion_std = np.std([mut_result[rep1][5], mut_result[rep2][5]])
            insert_deletion_yerr.append(insert_deletion_std)
        elif mut_result.__contains__(rep1):
            print("The group:",groupinfor.loc[idy].group, ": Rep2 is missing.")
            mutation.append(mut_result[rep1][1])
            mutation_yerr.append(0)
            insertO.append(mut_result[rep1][3])
            insertO_yerr.append(0)
            deletionO.append(mut_result[rep1][4])
            deletionO_yerr.append(0)
            insert_deletion.append(mut_result[rep1][5])
            insert_deletion_yerr.append(0)
        elif mut_result.__contains__(rep2):
            print("The group:", groupinfor.loc[idy].group, ": Rep1 is missing.")
            mutation.append(mut_result[rep2][1])
            mutation_yerr.append(0)
            insertO.append(mut_result[rep2][3])
            insertO_yerr.append(0)
            deletionO.append(mut_result[rep2][4])
            deletionO_yerr.append(0)
            insert_deletion.append(mut_result[rep2][5])
            insert_deletion_yerr.append(0)
        else:
            print("All repetitions in group:", groupinfor.loc[idy].group, " is missing.")
            mutation.append(0)
            mutation_yerr.append(0)
            insertO.append(0)
            insertO_yerr.append(0)
            deletionO.append(0)
            deletionO_yerr.append(0)
            insert_deletion.append(0)
            insert_deletion_yerr.append(0)

        if ck=='UNKNOWN':
            print("The group:",groupinfor.loc[idy].group, ": CK is missing.")
            ck_mutation.append(0)
            ck_insertO.append(0)
            ck_deletionO.append(0)
            ck_insert_deletion.append(0)
        else:
            ck_mutation.append(mut_result[ck][1])
            ck_insertO.append(mut_result[ck][3])
            ck_deletionO.append(mut_result[ck][4])
            ck_insert_deletion.append(mut_result[ck][5])


        glist.append(groupinfor.loc[idy].group)
        ck_glist.append(groupinfor.loc[idy].control)
    ## prepare for display

    ## print out pdf
    mutfile = os.path.join(output, 'mut_result.pdf')
    fig, (ax0, ax1) = plt.subplots(ncols=2, sharey=True)
    fig.set_size_inches(20, 9)
    width = 0.15
    bar1 = ax0.bar(groupinfor.index, mutation, width, color="#CC79A7")
    #bar1 = ax0.bar(groupinfor.index, replace, width, color='pink', yerr=replace_yerr, elinewidth=0.1, capsize=1.5)
    ax0.errorbar(groupinfor.index, mutation, yerr=mutation_yerr, fmt='', elinewidth=0.5, capsize=2, capthick=0.5, ls='None', ecolor='black')
    bar2 = ax0.bar(groupinfor.index + width, deletionO, width, color="#D55E00")
    #bar2 = ax0.bar(groupinfor.index + width, insertO, width, color='green', yerr=insertO_yerr, linewidth=0.5,capsize=1.5)
    ax0.errorbar(groupinfor.index+ width, deletionO, yerr=deletionO_yerr, fmt='', elinewidth=0.5, capsize=2, capthick=0.5, ls='None', ecolor='black')
    #bar3 = ax0.bar(groupinfor.index + width * 2, deletionO, width, color='blue', yerr=deletionO_yerr, linewidth=0.5,capsize=1.5)
    bar3 = ax0.bar(groupinfor.index + width * 2, insertO, width, color="#0072B2")
    ax0.errorbar(groupinfor.index + width * 2, insertO, yerr=insertO_yerr, fmt='', elinewidth=0.5, capsize=2, capthick=0.5,ls='None', ecolor='black')
    bar4 = ax0.bar(groupinfor.index + width * 3, insert_deletion, width, color="#009E73")
    #bar4 = ax0.bar(groupinfor.index + width * 3, insert_deletion, width, color='orange', yerr=insert_deletion_yerr,linewidth=0.5, capsize=1.5)
    ax0.errorbar(groupinfor.index + width * 3, insert_deletion, yerr=insert_deletion_yerr, fmt='', elinewidth=0.5, capsize=2,capthick=0.5,ls='None', ecolor='black')

    # ax.bar(reg.index, reg.delrate, color='blue')
    ax0.set_title('Treatment', size=15,fontdict = {'family': 'Times New Roman'})
    ax0.set_ylabel('All Mutation (%)', size=15,fontdict = {'family': 'Times New Roman'})
    ax0.set_xticks(groupinfor.index + 1.5 * width)
    #ax0.set_xticklabels(glist, rotation=35, size=6)
    ax0.set_xticklabels(glist, rotation=35, fontdict = {'family': 'Arial'}, size = 5)
    ax0.legend((bar1[0], bar2[0], bar3[0], bar4[0]), ('mutation_all', 'deletion_only', 'insert_only','insert&&deletion'))

    bar5 = ax1.bar(groupinfor.index, ck_mutation, width, color="#CC79A7")
    bar6 = ax1.bar(groupinfor.index + width, ck_deletionO, width, color="#D55E00")
    bar7 = ax1.bar(groupinfor.index + width * 2, ck_insertO, width, color="#0072B2")
    bar8 = ax1.bar(groupinfor.index + width * 3, ck_insert_deletion, width, color="#009E73")
    # ax.bar(reg.index, reg.delrate, color='blue')
    ax1.set_title('Control', size=15,fontdict = {'family': 'Times New Roman'})
    ax1.set_ylabel('All Mutation (%)', size=15,fontdict = {'family': 'Times New Roman'})
    ax1.set_xticks(groupinfor.index + 1.5 * width)
    #ax1.set_xticklabels(ck_glist, rotation=35, size=6)
    ax1.set_xticklabels(ck_glist, rotation=35, fontdict = {'family': 'Arial'}, size = 5)
    ax1.legend((bar5[0], bar6[0], bar7[0], bar8[0]), ('mutation_all', 'deletion_only', 'insert_only', 'insert&&deletion'))
    # plt.show()
    plt.savefig(mutfile, dpi=300, format="pdf")
    plt.close(fig)
    ## print out pdf
