import argparse
import sys
from CMlib import bwa, bwa_run, mut_rate, minimap2, minimap2_run
import os
import os.path
from pyfasta import Fasta
import pandas as pd
from subprocess import Popen, PIPE
import re

def main():
    args = check_options(get_options())
    build = args.build
    print("build:", build, "threads:", args.threads)

    # Build Gene_fasta index
    # This will generate a .fa.flat and .fa.gdx file in the same directory
    # as your input FASTA file

    print("## Step 1: Build index")
    indexfile = os.path.basename(args.genome)
    index = os.path.join(args.saved, indexfile+'.sa')
    if args.aligner == 'bwa':
        if build:
            bwa.index(args.bwa, args.genome, args.saved)
            print("bwa index build finished")
        else:
            print("Index already provided. Use", index)

        # Run bwa alignment
        print("## Step 2: Align to reference")
        print("Loading fastq files...")
        bwa_run.prepare(args.input, args.genome, args.saved, args.bwa, args.samtools, args.picard, args.threads)
        print("bwa mem finished!")
    elif args.aligner == 'minimap2':
        if build:
            print(args.minimap2)
            minimap2.index(args.minimap2, args.genome, args.saved)
            print("minimap2 index build finished")
        else:
            print("Index already provided. Use", index)

        print("## Step 2: Align to reference")
        print("Loading fastq files...")
        minimap2_run.prepare(args.input, args.genome, args.saved, args.minimap2, args.samtools, args.picard, args.threads)
        print("minimap2 alignment finished!")


    # Mutation rate calculation
    print("## Step 3: Calculate mutation rates")
    mut_rate.rate_cal(args.input, args.genome, args.result, args.saved)
    print("Mutation calculation finished!")



def check_options(parser):
    args = parser.parse_args()

    # Start check flash
    if args.flash:
        if not os.path.exists(args.flash):
            print("Can not locate flash, please input full path of flash\n")
            parser.print_help()
            sys.exit(1)
    else:
        flashpath = which('flash')
        if flashpath:
            flashversion=flash('flash')
            if flashversion == 'None':
                print("Can not locate flash, please input full path of flash\n")
                parser.print_help()
                sys.exit(1)
            else:
                args.flash = flashpath[0]
        else:
            print("Can not locate flash, please input full path of flash\n")
            parser.print_help()
            sys.exit(1)


    # Start check samtools
    if args.samtools:
        if not os.path.exists(args.samtools):
            print("Can not locate samtools, please input full path of samtools\n")
            parser.print_help()
            sys.exit(1)
    else:
        samtoolspath = which('samtools')
        if samtoolspath:
            samtoolsversion=samtools('samtools')
            if samtoolsversion == 'None':
                print("Can not locate samtools, please input full path of samtools\n")
                parser.print_help()
                sys.exit(1)
            else:
                args.samtools = samtoolspath[0]
        else:
            print("Can not locate samtools, please input full path of samtools\n")
            parser.print_help()
            sys.exit(1)
    # End check samtools

    # Start check picard
    if args.picard:
        if not os.path.exists(args.picard):
            print("Can not locate picard, please input full path of picard\n")
            parser.print_help()
            sys.exit(1)
    else:
        picardpath = which('picard')
        if picardpath:
            picardversion=picard('picard')
            if picardversion == 'None':
                print("Can not locate picard, please input full path of picard\n")
                parser.print_help()
                sys.exit(1)
            else:
                args.picard = picardpath[0]
        else:
            print("Can not locate picard, please input full path of picard\n")
            parser.print_help()
            sys.exit(1)
    # End check picard

    # Start check bwa
    if args.bwa:
        if not os.path.exists(args.bwa):
            print("Can not locate bwa, please input full path of bwa\n")
            parser.print_help()
            sys.exit(1)
        bwaversion = bwa.version(args.bwa)
        if bwaversion == 'None':
            print("Can not locate bwa, please input full path of bwa\n")
            parser.print_help()
            sys.exit(1)
    else:
        bwapath = which('bwa')
        if bwapath:
            bwaversion = bwa.version(bwapath[0])

            if bwaversion == 'None':
                print("Can not locate bwa, please input full path of bwa\n")
                parser.print_help()
                sys.exit(1)
            else:
                args.bwa = bwapath[0]
        else:
            print("Can not locate bwa, please input full path of bwa\n")
            parser.print_help()
            sys.exit(1)
    # End check bwa

    # Start check minimap
    if args.minimap2:
        if not os.path.exists(args.minimap2):
            print("Can not locate minimap2, please input full path of minimap2\n")
            parser.print_help()
            sys.exit(1)

        minimapversion = minimap2.minimap2version(args.minimap2)
        if minimapversion == 'None':
            print("Can not locate minimap, please input full path of minimap\n")
            parser.print_help()
            sys.exit(1)
    else:
        minimappath = which('minimap2')
        if minimappath:
            minimapversion = minimap2.version(minimappath[0])
            if minimapversion == 'None':
                print("Can not locate minimap, please input full path of minimap\n")
                parser.print_help()
                sys.exit(1)
            else:
                args.minimap2 = minimappath[0]
        else:
            print("Can not locate minimap, please input full path of minimap\n")
            parser.print_help()
            sys.exit(1)
    # End check minimap

    # Start check genome file
    if not os.path.exists(args.genome):
        print("Can not locate genome file, please input genome file.\n")
        parser.print_help()
        sys.exit(1)
    # End check genome file

    # Start check saved folder
    if not os.path.exists(args.saved):
        os.mkdir(args.saved)
    #End check saved folder

    # Start check result folder
    if not os.path.exists(args.result):
        os.mkdir(args.result)
    # End check result folder

    # Print Checked information
    print("#"*40)
    print("flash version:", args.flash, flashversion)
    print("bwa version:", args.bwa, bwaversion)
    print("minimap2 version:", args.minimap2, minimapversion)
    print("samtools version:", args.samtools, samtoolsversion)
    print("picard version:", args.picard, picardversion)
    print("genome file:", args.genome)
    print("input file:", args.input)
    print("tmp output folder:",  os.path.realpath(args.saved))
    print("result output folder:", os.path.realpath(args.result))
    print("threads number:", args.threads)
    print("#"*40)
    
    return args

def getch():
    """
    For yes/no choice
    """
    import sys, tty, termios
    fd = sys.stdin.fileno()
    old_settings = termios.tcgetattr(fd)
    try:
        tty.setraw(sys.stdin.fileno())
        ch = sys.stdin.read(1)
    finally:
        termios.tcsetattr(fd, termios.TCSADRAIN, old_settings)
    return ch

# Helper function for getting versions
def which(filename):
    """docstring for which"""
    locations = os.environ.get("PATH").split(os.pathsep)
    candidates = []
    for location in locations:
        candidate = os.path.join(location, filename)
        if os.path.isfile(candidate):
            candidates.append(candidate)
    return candidates


# Get versions
def samtools(filename):
    """
    :param filename:
    :return: samtools version
    """
    samtoolspath=which(filename)
    samtoolscmd = ' '.join([samtoolspath[0], '--version'])
    samtoolsrun = Popen(samtoolscmd, stdout=PIPE, stderr=PIPE, shell=True)
    i=samtoolsrun.stdout.readlines()[0]
    version = i.decode('utf-8').rstrip('\n')
    samtoolsrun.communicate()
    return version

def flash(filename):
    """
    :param filename:
    :return: flash version
    """
    flashpath=which(filename)
    flashcmd = ' '.join([flashpath[0], '--version'])
    #location= samtoolspath[0]
    flashrun = Popen(flashcmd, stdout=PIPE, stderr=PIPE, shell=True)
    i=flashrun.stdout.readlines()[0]
    version = i.decode('utf-8').rstrip('\n')
    flashrun.communicate()
    return version

def picard(filename):
    """
    :param filename:
    :return:
    """
    picardpath=which(filename)
    picardcmd = ' '.join([picardpath[0], 'ViewSam', '--version'])
    version = 'None'
    picardrun = Popen(picardcmd, stdout=PIPE, stderr=PIPE, shell=True)
    #print(picardcmd)
    for i in picardrun.stderr.readlines():

        i = i.decode('utf-8').rstrip('\n')

        if re.search('Version', i):
            (_, version) = i.split(':')
            print(version)

    picardrun.communicate()
    return version


def get_options():
    parser = argparse.ArgumentParser(description="CRISPRMatch is for location finding", prog='CRISPRMatch')
    parser.add_argument('--version', action='version', version='%(prog)s 1.0')
    parser.add_argument('-a', '--aligner', dest='aligner', help='preferred alignment tool', choices=['bwa','minimap2'], default='minimap2')
    parser.add_argument('-b', '--bwa', dest='bwa', help='bwa path')
    parser.add_argument('-m', '--minimap2', dest='minimap2', help='minimap2 path')
    parser.add_argument('-bd', '--build', dest='build', help='Build index or use provided?', choices=[True,False], default=True)
    parser.add_argument('-sm', '--samtools', dest='samtools', help='samtools path')
    parser.add_argument('-pi', '--picard', dest='picard', help='picard path')
    parser.add_argument('-f', '--flash', dest='flash', help='flash path')
    parser.add_argument('-g', '--genome', dest='genome', help='fasta format genome file', required=True)
    parser.add_argument('-i', '--input', dest='input', help='sample information input file', required=True)
    parser.add_argument('-s', '--save', dest='saved', help='tmp saved folder', default='tmpfiles')
    parser.add_argument('-r', '--result', dest='result', help='result saved folder', default='result')
    parser.add_argument('-t', '--threads', dest='threads', help='threads number or how may cpu you want to use',
                        default=1, type=int)
    parser.add_argument('--docker', default=False)

    return parser

if __name__ == "__main__":

    try:

        main()

    except KeyboardInterrupt:

        sys.stderr.write("User interrupt\n")

        sys.exit(0)
