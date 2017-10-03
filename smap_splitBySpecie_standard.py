#!/usr/local/bin/python

# Author : Remy Nicolle
# Date   : created 13 nov. 2015

# This script aims at getting a Graft specific BAM (usually Human) and a Host specific BAM (usually mice).


import sys
import pysam
import re
import getopt
import os
import subprocess
from datetime import datetime
import multiprocessing
import collections
import tempfile
import string
from copy import copy
from functools import partial


global TOTALFRAGMENTS
TOTALFRAGMENTS=0
global COUNTMULTIALIGNED
COUNTMULTIALIGNED=0

global RESCUESEQUENCE
RESCUESEQUENCE=0

global MORETHANONETRUEPAIR
MORETHANONETRUEPAIR=0

global MULTIALIGNEDAMBIGUITY
MULTIALIGNEDAMBIGUITY=0

global TRUEPAIRNOTPRIMARYALIGN
TRUEPAIRNOTPRIMARYALIGN=0

global TRUEPAIRSECONDARY
TRUEPAIRSECONDARY =0

global NOTALIGNED
global GRAFTONLY
global HOSTONLY
global WRONGPAIR
global AMBIGUOUS
global MULTIMAPPEDNOTDEALTWITH

NOTALIGNED="notAligned"
GRAFTONLY="graftonly"
HOSTONLY="hostonly"
WRONGPAIR="wrongpair"
AMBIGUOUS="ambiguous"
MULTIMAPPEDNOTDEALTWITH="multimappedNotDealtWith"
MULTIMAPPEDHOST="multimappedhost"
MULTIMAPPEDGRAFT="multimappedgraft"

#---------------------------------------------------------------------------------#
# create a BAM writer for a given specie and based on a multispecie  input bam
def prepareBAMwriter(originalBAM,outputfilepath,specieIDpattern,rehead):
    olddicheader=originalBAM.header
    speciedicheader=copy(originalBAM.header)
    if rehead :
        speciedicheader["SQ"]=[]
        for i in range(0,len(olddicheader["SQ"])):
            if(  olddicheader["SQ"][i]["SN"].startswith(specieIDpattern)):
                toaddseq=copy(olddicheader["SQ"][i])
                toaddseq["SN"]=toaddseq["SN"].replace(specieIDpattern,"")
                speciedicheader["SQ"].append(toaddseq)
    return pysam.AlignmentFile(outputfilepath, mode="wb", header=speciedicheader)

#---------------------------------------------------------------------------------#
# stripping out name of sample from file by removing the last .*
def getSampleFileName(origfile):
    samplist=origfile.split(".")
    samplist.pop()
    sampnam=''
    for num in samplist:
        sampnam +=num +"."
    return sampnam[:-1]







#-----------------------------------------------------------------------------------#
# bamwriter new reference
def bamwriter2referenceseqs(abamwriter) :
    newref=[]
    for sq in copy(abamwriter.header)["SQ"] :
        newref.append(sq["SN"])
    return newref

#-----------------------------------------------------------------------------------#
# write a sequence in a bam with a new set of reference sequences
def writeNewRef(listofreads,awriter,anewreflist,specieIDpattern):
    global RESCUESEQUENCE
    for r in listofreads :
        if ((r.seq==None)|(r.qual==None)):
                RESCUESEQUENCE+=1
        newname=r.reference_name.replace(specieIDpattern,"")
        newmatename=r.next_reference_name.replace(specieIDpattern,"")
        r.rname=anewreflist.index(newname)
        r.rnext=anewreflist.index(newmatename)
        awriter.write(r)



#---------------------------------------------------------------------------------#
#Verigyinf file exists and is a bam file
def verifyBam(afilename,what):
    if afilename=="":
        print afilename +" file, serving as "+ what+", does not exist"
        sys.exit()
    try:
        tmp=pysam.AlignmentFile(afilename,"rb")
    except Exception as e:
        print afilename +" file, serving as "+ what+", is not a BAM file"
        print "error :"
        print e
        sys.exit()
    if(not tmp.is_bam):
        print afilename +" file, serving as "+ what+", is not a BAM file"
        sys.exit()
    tmp.close()
    del tmp

#---------------------------------------------------------------------------------#
#Keeping track of time
def keepingTrack(what,info,sofar,verbose):
    sofar+= datetime.now().strftime("%d %a %b (20%y) %H:%M:%S")+ "    "+what+"    "+info+"\n"
    if verbose:
        print info
    return sofar



#---------------------------------------------------------------------------------#
# test wether one read considers the other as its pair
def areAsymetricPairs(read1,read2):
    if(read1.is_unmapped |read2.is_unmapped):
        return False
    else:
        return(((read1.reference_name == read2.next_reference_name) &
        (read1.reference_start == read2.next_reference_start))|
        ((read1.next_reference_name == read2.reference_name) &
        (read1.next_reference_start == read2.reference_start ))&
        (read1.qname==read2.qname)&
        ((read1.is_read1 & read2.is_read2)|(read1.is_read2 & read2.is_read1))&
        (read1.is_paired) & (read2.is_paired)
        )

def areTruePairs(read1,read2):
    if(read1.is_unmapped |read2.is_unmapped):
        return False
        #    elif( (not read1.is_proper_pair )|(not read2.is_proper_pair)):
        #        return False
    elif ( (not read1.is_paired) | (not read2.is_paired)):
        return False
    else:
        return(((read1.reference_name == read2.next_reference_name) &
                (read1.reference_start == read2.next_reference_start))&
               ((read1.next_reference_name == read2.reference_name) &
                (read1.next_reference_start == read2.reference_start ))&
               (read1.qname==read2.qname)&
               ((read1.is_read1 & read2.is_read2)|(read1.is_read2 & read2.is_read1))
               )


#---------------------------------------------------------------------------------#
# Processing pairs of reads for two species and with dealing with multiple alignement
#
#   Specifc behavior of this function :
#   - This function will only rely on the Alignement Score given as an extra BAM field (usually AS)
#   - This function will look at all asymetric pairs (at least one read considers the other as its mate)
#   - only true mates will be found in the output bams (properely paired AND symetric pairs were both consider the other as its mate)
#   - only non ambiguous true mates (mapped on same organism) that have higher alignment scores than any other asymetric pairs wich might map to another organism are selected

# subalignementScoreField = "XS" and alignementScoreField="AS"
def processBundle2speciesPairEndMultiAlignedWithScores(bundle,graftpattern,hostpattern,graftWriter,hostWriter,graftreferences,hostreferences,alignementScoreField,outputmultipleAlignments):


    global MORETHANONETRUEPAIR
    global MULTIALIGNEDAMBIGUITY
    
    global TRUEPAIRNOTPRIMARYALIGN
    global TRUEPAIRSECONDARY

    if(len(bundle)>2):
        bundle=list(set(bundle))
    
    

    # saving sequence in case its not stored in duplicare alignments
    read1Seq =''
    read2Seq =''
    read1Qual=''
    read2Qual=''
    for b in bundle:
        if(b.is_read1):
            if(b.seq!=None):
                if(len(b.seq)>len(read1Seq)):
                    read1Seq=b.seq
            if(b.qual!=None):
                if(len(b.qual)>len(read1Qual)):
                    read1Qual=b.qual
        if(b.is_read2):
            if(b.seq!=None):
                if(len(b.seq)>len(read2Seq)):
                    read2Seq=b.seq
            if(b.qual!=None):
                if(len(b.qual)>len(read2Qual)):
                    read2Qual=b.qual


    # bundle of sequences to bundle of pairs !
    # stores all possible asymetric pairs (at least one read considers the other as its pair) in "allpossiblepairedbundle"
    # and all true pairs (both consider each other as their mate) in "truepairs"
    allpossiblepairedbundle=[]
    truepairs=[]
    for b in bundle :
        if( b.is_read1 & b.is_paired):
            for c in bundle :
                if c.is_read2 :
                    ispair=False
                    ispair=areAsymetricPairs(b,c)
                    if areTruePairs(b,c):
                        truepairs.append([b,c])
                    if( ispair):
                        allpossiblepairedbundle.append([b,c])

# if these are empty, just stop here, the reads are not aligned correctly or not at all
    if(len(allpossiblepairedbundle)==0):
        return NOTALIGNED
    if(len(truepairs)==0):
        return NOTALIGNED


    # now, we get the maximum alignment score of all asymetric pairs.
    # And we store the asymetric pairs that have this alignment score
    # the alignment score of a pair is here the sum of the alignement score of each read
    # getting all asymetric pairs with max score
    maxasymetricscore=float("-inf")
    bestasymetricpairedbundles=[]
    for p in allpossiblepairedbundle :
        thiscore=p[0].get_tag(alignementScoreField)+p[1].get_tag(alignementScoreField)
        if thiscore>maxasymetricscore:
            bestasymetricpairedbundles=[]
            bestasymetricpairedbundles.append(p)
            maxasymetricscore=thiscore
        elif thiscore==maxasymetricscore:
             bestasymetricpairedbundles.append(p)

    # and we do the same for true pairs but only if the maximum true pair score is at last as big as the maxasymetric pair score.
    maxtruepairscore=maxasymetricscore
    besttruepairs=[]
    for p in truepairs :
        thiscore=p[0].get_tag(alignementScoreField)+p[1].get_tag(alignementScoreField)
        if thiscore>maxtruepairscore:
            besttruepairs=[]
            besttruepairs.append(p)
            maxtruepairscore=thiscore
        elif thiscore==maxtruepairscore:
            besttruepairs.append(p)

    # This is a security check ... If this is empty, we probably have a problem somewhere ...
    # or no true pair has a higher score than an asymetric pair, so weve got nothing to write
    if( (len(bestasymetricpairedbundles)==0 )| (len(besttruepairs)==0)):
        return WRONGPAIR


#Now, we've got true pairs (the ones we could keep and write into a result BAM) that have the maximum possible score for this read pairs.
# And we also have the asymetric pairs with their maximum possible scores, these contain the best possible true pairs
# if you are here, then the true pairs have a score equal or higher than the best asymetrical pairs


    allbestreads =set()
    for b in bestasymetricpairedbundles :
        allbestreads.add(b[0])
        allbestreads.add(b[1])

# first, check if you have best pairs (asymetrical or true) that are ambiguous
    # for these BEST asymetric alignment pairs, if we have both graft and host, then this is an ambiguous mate
    # reminder : the true pairs are also in the asymetric pairs
    hasgraft=False
    hashost=False
    for b in allbestreads:
        hasgraft =hasgraft| b.reference_name.startswith(graftpattern)
        hashost =hashost| b.reference_name.startswith(hostpattern)

    if(hashost & hasgraft):
        return AMBIGUOUS

    if((len(bestasymetricpairedbundles)>1)&(hashost & hasgraft)):
        MULTIALIGNEDAMBIGUITY+=1

    if(len(besttruepairs)>1):
        MORETHANONETRUEPAIR+=0
    if((besttruepairs[0][0].is_secondary)|(besttruepairs[0][1].is_secondary)):
       TRUEPAIRNOTPRIMARYALIGN+=1
    if((besttruepairs[0][0].is_supplementary)|(besttruepairs[0][1].is_supplementary)):
          TRUEPAIRSECONDARY+=1
       
    readstowrite=set()
    if outputmultipleAlignments :
        readstowrite=allbestreads
    elif len(besttruepairs)==1 :
        readstowrite=set(besttruepairs[0])
    else :
        readstowrite=set()
        for b in besttruepairs:
            readstowrite.add(b[0])
            readstowrite.add(b[1])



    for b in readstowrite:
        if((b.seq==None)|(b.qual==None)):
            if(b.is_read1):
                b.seq=read1Seq
                b.qual=read1Qual
            elif(b.is_read2):
                b.seq=read2Seq
                b.qual=read2Qual


    if(hashost & hasgraft):
        return AMBIGUOUS
    elif hashost:
        writeNewRef(readstowrite,hostWriter,hostreferences,hostpattern)
        return HOSTONLY
    elif hasgraft:
        writeNewRef(readstowrite,graftWriter,graftreferences,graftpattern)
        return GRAFTONLY
    else:
        print "WRONG WRONG WRONG SOMETHING WENT WRONG WRONG WRONG totoudidou ..."
        sys.exit(1)




#---------------------------------------------------------------------------------#
# Processing pairs of reads for two species and NOT DEALING with multiple alignement
def processBundle2speciesPairEndWithoutMultiAligned(bundle,graftpattern,hostpattern,graftWriter,hostWriter,alignementScoreField,subalignementScoreField):
    if(len(bundle)==2):# if two reads, probably paired reads
    #all this is with a bundle of two reads
        if((not  bundle[0].aligned )|( not bundle[1].aligned )):# if one is not aligned, the pair is not aligned
            return notAligned
        elif(not (((bundle[0].pe_which=="first") &( bundle[1].pe_which=="second"))|
             ((bundle[1].pe_which=="first") & (bundle[0].pe_which=="second")))):# if got two reads but not mates
             return WRONGPAIR
        elif(bundle[0].iv.chrom.startswith(graftpattern) & bundle[1].iv.chrom.startswith(graftpattern)): # if mates are grafts
        #from here its is pure GRAFT
            if((bundle[0].optional_field(alignementScoreField) > bundle[0].optional_field(subalignementScoreField))
               |(bundle[1].optional_field(alignementScoreField) > bundle[1].optional_field(subalignementScoreField))):# if this score is better than subotpiomal score
                    if( (bundle[0].proper_pair )&( bundle[1].proper_pair)):#if not properly paired
                        for a in bundle:
                            graftWriter.write( a) # GREAT, this is graft write to graft bam
                        return GRAFTONLY
                    else:# then this is a wrong pair
                        return WRONGPAIR
            else:#then this is an ambiguous read (suboptimal score is the same)
                return MULTIMAPPEDNOTDEALTWITH
        elif(bundle[0].iv.chrom.startswith(hostpattern) & bundle[1].iv.chrom.startswith(hostpattern)): #both mapped to host
            #from here its pure HOST
            if((bundle[0].optional_field(alignementScoreField) > bundle[0].optional_field(subalignementScoreField)) |
               (bundle[1].optional_field(alignementScoreField) > bundle[1].optional_field(subalignementScoreField))):# if this score is better than subotpiomal score
                if( (bundle[0].proper_pair )&( bundle[1].proper_pair)):
                    for a in bundle:
                        hostWriter.write( a)
                    return HOSTONLY
                else:
                    return WRONGPAIR
            else:
                return MULTIMAPPEDNOTDEALTWITH
        else:
            return AMBIGUOUS
    elif(len(bundle)<2):
        return WRONGPAIR
    else:
        return MULTIMAPPEDNOTDEALTWITH




#---------------------------------------------------------------------------------#
#  single thread of processing a sorted BAM

def htseq_bundle_multiple_alignments( sequence_of_alignments ):
    """Some alignment programs, e.g., Bowtie, can output multiple alignments,
    i.e., the same read is reported consecutively with different alignments.
    This function takes an iterator over alignments and bundles consecutive
    alignments regarding the same read to a list of Alignment objects and
    returns an iterator over these.
    """
    alignment_iter = iter( sequence_of_alignments )
    algnt = alignment_iter.next()
    ma = [ algnt ]
    for algnt in alignment_iter:
        if algnt.qname != ma[0].qname:
            yield ma
            ma = [ algnt ]
        else:
            ma.append( algnt )
    yield ma





def singleThreadBamSplit(inputbam,outputprefix,graftpattern,hostpattern,alignementScoreField,verbose,log,outputmultipleAlignments):
    global TOTALFRAGMENTS
    global COUNTMULTIALIGNED

    graftbamfilename =outputprefix+"_"+graftpattern+".bam"
    graftWriter=prepareBAMwriter(inputbam,graftbamfilename,graftpattern,True)
    graftreferences=bamwriter2referenceseqs(graftWriter)

    
    hostbamfilename =outputprefix+"_"+hostpattern+".bam"
    hostWriter=prepareBAMwriter(inputbam,hostbamfilename,hostpattern,True)
    hostreferences=bamwriter2referenceseqs(hostWriter)
    
    labels=[]


    for bundle in htseq_bundle_multiple_alignments( inputbam ):
        TOTALFRAGMENTS+=1
        if len(bundle)>2:
            COUNTMULTIALIGNED+=1
        if(( TOTALFRAGMENTS /1000000 ==  TOTALFRAGMENTS /1000000.0 )):
            log=keepingTrack("Step 2 : splitting",str(TOTALFRAGMENTS) +" fragments processed.",log,verbose)
            if(COUNTMULTIALIGNED==0):
                print "No multialigned fragments found so far ("+str(TOTALFRAGMENTS)+"processed). Check the mapping was ran in multalignment mode."

        labels.append(processBundle2speciesPairEndMultiAlignedWithScores(bundle,graftpattern,hostpattern,graftWriter,hostWriter,graftreferences,hostreferences,alignementScoreField,outputmultipleAlignments))
    
    hostWriter.close()
    graftWriter.close()
    inputbam.close()
    return collections.Counter(labels)



#---------------------------------------------------------------------------------#
#
def checkTag(inputfile,atag):
    try:
        tmp=pysam.AlignmentFile(inputfile,"rb")
        a=tmp.next().get_tag(atag)
        return atag
    except Exception as e:
        print "The given alignment score BAM field ("+atag+") does not exist in this bam file"
        sys.exit(0)

def checkPattern(inputfile,pattern,what,arg):
    H=pysam.AlignmentFile(inputfile,"rb").header["SQ"]
    found=False
    for h in H :
        found= h["SN"].startswith(pattern)
        if found:
            break
    if (not found):
        print "Pattern "+pattern+" used as the "+what+" (argument "+arg+"), was not found among reference sequences."
        sys.exit(0)



#---------------------------------------------------------------------------------#
#       MAIN FUNCTION : one function to run them all
#---------------------------------------------------------------------------------#

def main(debug):
    #Message to be displayed in case help is needed to run the exectuable
    helpstr="Usage : smapFilter.py  -i <inputBAMfile> [-g <graftSpecieID>] [-h <hostSpecieID>] [-o <outputBamFilesPrefix>] [-k -v] [-f <alignementScoreFieldInBAM>]\n\n"
    helpstr+=" -i inputBAMfile.bam  : A bam file coming from the simultaneous mapping of a FASTQ"
    helpstr+="   on two (or more?) reference genome (i.e. mouse and human for a xenograft sample experiment)\n"
    helpstr+=" -g hs  : a pattern corresponding to the beginning of chromosome name of the graft specie. (Default hs)\n"
    helpstr+=" -h mmu   : a pattern corresponding to the beginning of chromosome name of the host specie. (default mmu)\n"
    helpstr+=" optional arguments :\n"
    helpstr+=" -n <NumberOfThreads> : Number of threads to run in parallel on a single multi core machine. Default is 1.\n"
    helpstr+=" -k : keep the lexicographically sorted BAM, default is to remove it\n"
    helpstr+=" -f AlignementScoreField : The name of the alignement score field (default is AS)\n"
    helpstr+=" -o : output file prefix (including  absolute or conditional path if needed)\n"
    helpstr+=" -v : Verbose. Print information during process\n"
    helpstr+="  --help : displays this message"
    #general variables
    log=""
    # setting default parameters
    inputfile = "wholeSamDoubleMapBWAWithMultiHit_downsample.bam" # input bam file
    graftpattern ="hs"
    hostpattern =  "mmu"
    removesortedbam=True
    outputprefix=""
    alignementScoreField="AS"
    verbose=False
    nthreads=1
    outputmultipleAlignments=False
    
    print verbose
    #---------------------------------------------------------------------------------#
    # retrieving command line parameters
    if(not debug):
        verbose=False
        try:
            opts, args = getopt.getopt(sys.argv[1:], "i:g:h:kvf:m",["--help",]);
        except getopt.GetoptError:
            print helpstr
            sys.exit(2);
        if len(sys.argv)==1:
            print helpstr
            sys.exit()
        for opt, arg in opts:
            if opt =="--help":
                print helpstr
                sys.exit()
            elif opt =="-i":
                inputfile = arg
            elif opt =="-g":
                graftpattern = arg
            elif opt =="-h":
                hostpattern = arg
            elif opt =="-k":
                removesortedbam = False
            elif opt=="-f":
                alignementScoreField=arg
            elif opt=="-v":
                verbose=True
            elif opt=="-m":
                outputmultipleAlignments=True
            elif opt =="-o":
                outputprefix=arg
    
    #---------------------------------------------------------------------------------#
    #starting, input verification and file names process
    log=keepingTrack("start","",log,verbose)
    verifyBam(inputfile,"original BAM file to process")

    checkTag(inputfile,"AS")
    checkPattern(inputfile,graftpattern,"graftpattern","-g")
    checkPattern(inputfile,hostpattern,"hostpattern","-h")

    # dealing with sample names and files
    samplefilename=getSampleFileName(inputfile)
    sortedbamname=samplefilename+"_sorted"
    sortedbamfil=sortedbamname+".bam"

    if(outputprefix==""):
        outputprefix=samplefilename

    #---------------------------------------------------------------------------------#
    # step 1 sorting

    log=keepingTrack("Step 1 : sorting","Lexicographically sorting the input BAM",log,verbose)
    pysam.sort("-n", inputfile,sortedbamname)
    print sortedbamfil
    verifyBam(sortedbamfil,"Sorted BAM")
    log=keepingTrack("Step 1 : sorting","Sorting done",log,verbose)
    #---------------------------------------------------------------------------------#
    # step 2 spliting by specie
    log=keepingTrack("Step 2 : splitting","Start splitting per species",log,verbose)
    inputbam = pysam.AlignmentFile(sortedbamfil,"rb" )

    stats=singleThreadBamSplit(inputbam,outputprefix,graftpattern,hostpattern,alignementScoreField,verbose,log,outputmultipleAlignments)

    log=keepingTrack("Step 2 : splitting","Done splitting per species",log,verbose)
    print log
    #---------------------------------------------------------------------------------#
    # writing results and cleaning up
    global TOTALFRAGMENTS

    if(removesortedbam):
        os.remove(sortedbamfil)

    statfile=open(outputprefix+".stat.txt","w")
    statfile.write("totalReads"+"\t"+str(TOTALFRAGMENTS)+"\n")
    for k in stats.keys():
        statfile.write(k+"\t"+str(stats[k])+"\n")
    statfile.close()

            
    global MULTIALIGNEDAMBIGUITY
    print "Ambiguity found using multiple alignments:"
    print MULTIALIGNEDAMBIGUITY
            
    
    global RESCUESEQUENCE
    print "Missing sequences:"
    print  RESCUESEQUENCE

    
    global MORETHANONETRUEPAIR
    print "Number of times more than one true pair was present :"
    print MORETHANONETRUEPAIR
   
    global TRUEPAIRNOTPRIMARYALIGN
    print "Number of times more the one true pair was not the primary alignment :"
    print TRUEPAIRNOTPRIMARYALIGN
       
    global TRUEPAIRSECONDARY
    print "Number of times more the one true pair was a supplementary alignment :"
    print TRUEPAIRSECONDARY

    global COUNTMULTIALIGNED
    print "Number of multialigned read pairs :"
    print COUNTMULTIALIGNED


if __name__ == "__main__":
    main(False)


#    print datetime.now().strftime("%d %a %b (20%y) %H:%M:%S")+ "    info    "+ str(counter["totalpairs"]) + " reads"
#end_time = datetime.now()
#print('Total Duration: {}'.format(end_time - start_time))


#  inputfile = pysam.AlignmentFile('wholeSamDoubleMapBWAWithMultiHit_downsample.bam','rb')
