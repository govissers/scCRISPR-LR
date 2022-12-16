import os
from subprocess import Popen
from subprocess import PIPE
import re
import shutil
import pandas as pd

# Rewrite the FASTA file to contain the deletion variants
def generate(info, ref):
    # Read in geneinfo file
    sampleinfo = pd.read_table(info, index_col='Index')

    # Create a copy of the existing fasta file
    refnoext = os.path.splitext(ref)[0]
    outname = refnoext+'_lgdel.fa'

    # Write to a new file
    infile = open(ref, 'r')
    outio = open(outname, 'w')
    print('Generating deletion constructs for:')
    for idx in sampleinfo.index:
        genefile = os.path.dirname(os.path.abspath(info)) + '/' + sampleinfo.loc[idx]['Gene_File'] + '.txt'
        genetable = pd.read_table(genefile, index_col="geneid")
        lines = infile.readlines()
        for i in range(0, len(lines) - 1):
            line = lines[i]
            nextline = lines[i + 1]
            if line[0] == '>':
                gene = line[1:]
                gene = gene.replace("\n", '')
                print(gene)
                if i == len(lines) - 2:
                    nextline = nextline+"\n"
                outio.write(line)
                outio.write(nextline)
                outio.write('\n')

                # Get cutsites
                cut1pos = genetable.loc[gene]['start_1']
                cut2pos = genetable.loc[gene]['start_2']
                cut3pos = genetable.loc[gene]['start_3']
                
                if cut2pos is not None:
                    # Write in large deletions
                    outio.write('>'+gene+'_DEL1\n')
                    outio.write(nextline[0:cut1pos] + nextline[cut2pos:])
                    outio.write('\n')

                    if cut3pos is not None:
                        outio.write('>'+gene+'_DEL2\n')
                        outio.write(nextline[0:cut2pos] + nextline[cut3pos:])
                        outio.write('\n')

                        outio.write('>'+gene+'_DEL3\n')
                        outio.write(nextline[0:cut1pos] + nextline[cut3pos:])
                        outio.write('\n')
    return outname

                

