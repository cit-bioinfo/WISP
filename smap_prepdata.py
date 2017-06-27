#!/usr/local/bin/python

# Author : Remy Nicolle
# Date   :  created 8 jan. 2016

# This script aims at preparing genomic datasets (genome fasta files and transcriptome gtf files)
# for the SMAP analysis


#  python /datacit/Development/RNAseq_YR/dev-pdx-git/smap_prepdata.py  --graftfasta /datacit/00_DATABANKS/ensembl_hg19/Homo_sapiens.GRCh37.75.dna.primary_assembly.fa --hostfasta /datacit/00_DATABANKS/ensembl_hg19/Homo_sapiens.GRCh37.75.dna.primary_assembly.fa --graftgtf /datacit/00_DATABANKS/ensembl_hg19/gencode.v19.annotation.gtf --hostgtf /datacit/00_DATABANKS/ensembl_hg19/Homo_sapiens.GRCh37.75.gtf -o ~/

import sys
import re
import getopt
import os
import subprocess
from datetime import datetime
import collections
import string
from copy import copy
from functools import partial


#---------------------------------------------------------------------------------#
#Keeping track of time
def keepingTrack(what,info,sofar,verbose):
    sofar+= datetime.now().strftime("%d %a %b (20%y) %H:%M:%S")+ "    "+what+"    "+info+"\n"
    if verbose:
        print info
    return sofar

#---------------------------------------------------------------------------------#
#Testing wether the GTF has the pattern
def verifyGTF(gtfile,graftPattern,hostPattern,fileUse):
    hasHostPattern=False
    hasGraftPattern=False

    try:
       f = open(gtfile)
    except EnvironmentError:
        print "Could not access or read the "+fileUse+" file."
        sys.exit(0)
    
    with f:
        for lines in f:
            if not lines.startswith("#"):
                if lines.strip().split("\t")[0].find(hostPattern) != -1:
                    hasHostPattern=True
                    break
                if  lines.strip().split("\t")[0].find(graftPattern) != -1:
                    hasGraftPattern=True
                    break
    f.close()
    if hasHostPattern | hasGraftPattern :
        warnin="The pattern that should specifically identify chromosomes from"
        if (hasHostPattern) &(not hasGraftPattern):
            warnin+= " the host specie"
        elif (hasGraftPattern) &(not hasHostPattern):
            warnin+= " the graft specie"
        else :
            warnin+= "  both the host and the graft species"
        warnin+= "  were found in the original "+fileUse+" file.\n"
        warnin+= "To change the host pattern, add the argument -h (or -g) followed by a new simple pattern while running the script in command line (example : python smap_prepdata.py -h homosapiens)."
        print warnin
        sys.exit(0)


def verifyFasta(fastafile,graftPattern,hostPattern,fileUse):
    hasHostPattern=False
    hasGraftPattern=False
    try:
        f = open(fastafile)
    except EnvironmentError:
        print "Could not access or read the "+fileUse+" file."
        sys.exit(0)

    with f:
        for lines in f:
            if  lines.startswith(">"):
                if lines.strip().split("\t")[0].find(hostPattern) != -1:
                    hasHostPattern=True
                    break
                if  lines.strip().split("\t")[0].find(graftPattern) != -1:
                    hasGraftPattern=True
                    break
    f.close()
    if hasHostPattern | hasGraftPattern :
        warnin="The pattern that should specifically identify chromosomes from"
        if (hasHostPattern) &(not hasGraftPattern):
            warnin+= " the host specie"
        elif (hasGraftPattern) &(not hasHostPattern):
            warnin+= " the graft specie"
        else :
            warnin+= "  both the host and the graft species"
        warnin+= "  were found in the original "+fileUse+" file.\n"
        warnin+= "To change the host pattern, add the argument -h (or -g) followed by a new simple pattern while running the script in command line (example : python smap_prepdata.py -h homosapiens)."
        print warnin
        sys.exit(0)



def main():
    #Message to be displayed in case help is needed to run the exectuable
    helpstr="Usage : smap_prepdata.py  --graftfasta <graftFasta.fa> --hostfasta <hostFasta.fa> --graftgtf <graftGTF.gtf> --hostgtf <hostGTF.gtf> -o <outputDirectory> [-s -b -g <graftSpecificPattern> -h <hostSpecificPatternn> ]\n\n"

    helpstr+=" --graftfasta graftFasta.fa  : fasta file of the graft organism\n"
    helpstr+=" --hostfasta hostFasta.fa   : fasta file of the host organism\n"
    helpstr+=" --graftgtf graftGTF.gtf   : gtf file of the graft organism\n"
    helpstr+=" --hostgtf hostGTF.gtf    : gtf file of the host organism\n"
    helpstr+=" -o graftFasta.fa  : Directory to store all"
    helpstr+=" optional arguments \n"
    helpstr+=" -s                           : Create a STAR index, necessary to run a STAR alignment (Will not run by default)\n"
    helpstr+=" -b                           : Create a BWA index, necessary to run a BWA alignment (Will not run by default)\n"
    helpstr+=" -g <graftSpecificPattern>    : fasta file of the graft organism (hs by default)\n"
    helpstr+=" -h <hostSpecificPattern>     : fasta file of the host organism (mmu by default)\n"
    helpstr+=" -v                           : run in verbose mode to display a little information\n"
    helpstr+="  -h                      : displays this message\n"

    #general variables
    log=""
    # setting default parameters
    graftFastaPath=""
    hostFastaPath=""
    graftGTFPath=""
    hostGTFPath=""
    outputDirectory=""
    verbose=False
    graftPattern="hs"
    hostPattern="mmu"
    doBWAIndex=False
    doStarIndex=False
    verbose=False
    

    #---------------------------------------------------------------------------------#
    # retrieving command line parameters

    try:
        opts, args = getopt.getopt(sys.argv[1:], "o:sbh",["graftfasta","hostfasta","graftgtf","hostgtf"]);
        print opts
        print args

    except getopt.GetoptError:
        print "opterror"
        print helpstr
        sys.exit(2);
    if len(sys.argv)==1:
        print helpstr
        sys.exit()
    for opt, arg in opts:
        if opt =="-h":
            print helpstr
            sys.exit(0)
        elif opt =="--graftfasta":
            graftFastaPath = arg
        elif opt =="--hostfasta":
            hostFastaPath = arg
        elif opt =="--graftgtf":
            graftGTFPath = arg
        elif opt =="--hostgtf":
            hostGTFPath = False
        elif opt=="-o":
            outputDirectory=arg
        elif opt=="-s":
            doStarIndex=True
        elif opt=="-b":
            doBWAIndex=True
        elif opt=="-h":
            graftPattern=True
        elif opt=="-g":
            hostPattern=True
        elif opt=="-v":
            verbose=True

    if graftFastaPath=="":
        print "Missing graft Fasta file"
        sys.exit()
    if hostFastaPath=="":
        print "Missing host Fasta file"
        sys.exit()
    if graftGTFPath=="":
        print "Missing graft GTF file"
        sys.exit()
    if hostGTFPath=="":
        print "Missing host GTF file"
        sys.exit()
    if outputDirectory=="":
        print "Missing output directory"
        sys.exit()


    #---------------------------------------------------------------------------------#
    #starting, input verification and file names process
    log=keepingTrack("start","",log,verbose)
    log=keepingTrack("Checking input files","host gtf file",log,verbose)
    verifyGTF(hostGTFPath,graftPattern,hostPattern,"host GTF")
    log=keepingTrack("Checking input files","graft gtf file",log,verbose)
    verifyGTF(graftGTFPath,graftPattern,hostPattern,"graft GTF")
    log=keepingTrack("Checking input files","host fasta file",log,verbose)
    verifyFasta(hostFastaPath,graftPattern,hostPattern,"host Fasta")
    log=keepingTrack("Checking input files","graft fasta file",log,verbose)
    verifyFasta(graftFastaPath,graftPattern,hostPattern,"graft Fasta")


if __name__ == "__main__":
    main()