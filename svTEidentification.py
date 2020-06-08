#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
	Identifying complete TE differences between two assemblies using Assemblytics/RaGOO output
	==========================
	:author:  François Sabot
	:contact: francois.sabot@ird.fr
	:date: 10/10/2019
	:version: 0.1
	:licence: GPLv3
	Script description
	------------------
	svTEidentification.py will inform about the new TE putative insertion based on a Assemblytics BED file.
    A TE database (mutlifasta file) must be provided
    Requires bedtools and blast to be accessibl in the path
	-------
	>>> svTEidentification.py -i Assemblytics.bed -d TEdatabase -o output

	Help Programm
	-------------
	information arguments:
		- \-h, --help
						show this help message and exit
		- \-v, --version
						display svTEidentification.py version number and exit
	Input mandatory infos for running:
		- \-i <filename>, --input <filename>
						BED file issued from assemblytics analysis
		- \-r <filename>, --reference <filename>
						reference sequence file used in assemblytics
        - \-a <filename>, --alternate <filename>
						alternate sequence file used in assemblytics
        - \- d <filename>, --database <filename>
						TE database in multifasta format
		- \-o <filename>, --out <filename>
						Prefix of output files
		- \-s <minimalPercentage>, --size <minimalPercentage>
						Minimal percentage of identity and size for a TE hit to be conserved (Optional, default 90)
"""

#Import
import os
current_dir = os.path.dirname(os.path.abspath(__file__))+"/"

import subprocess as sp

## Python modules
import argparse
from Bio import SeqIO
#import numpy as np
#import pandas as pd
#import matplotlib.pyplot as plt

##################################################
## Variables Globales
##################################################
version="0.1"
VERSION_DATE='10/10/2019'
debug="False"
#debug="True"

##################################################
## Functions
##################################################

def relativeToAbsolutePath(relative):
	if relative[0] != "/":			# The relative path is a relative path, ie do not starts with /
		command = "readlink -m "+relative
		absolutePath = sp.check_output(command, shell=True).decode("utf-8").rstrip()
		return absolutePath
	else:							# Relative is in fact an absolute path, send a warning
		absolutePath = relative;
		return absolutePath

def launchingShell(command): # Function to launch and display command
    '''Launcher for command'''
    #print(command)
    process = sp.Popen(command, stdout=sp.PIPE, shell=True, stderr=sp.PIPE)
    (outputCommand, errorCommand) = process.communicate()
    return (outputCommand, errorCommand)


def indelExtraction(inputFile,outputPrefix):
    '''Transform a global Assemblytics BED in two BED, one with insertion, one with deletion'''
    print ("Extraction of informations for Deletions & Insertions")
    grepCommand = "cat " + inputFile + "  | grep Deletion > " + outputPrefix + "_deletion.bed && cat " + inputFile + "  | grep Insertion > " + outputPrefix + "_insertion.bed"
    launchingShell(grepCommand)

def BEDconverter(inputBed,outputBed):
    '''Transform an Assemblytics BED in a BED based on the alternate sequence'''
    print ("Converting BED from reference to alternate")
    outputHandle = open(outputBed, "w")
    with open(inputBed, "r") as inputHandle:
        for line in inputHandle:
            if line.startswith("reference"):#Passing through header
                continue
            mainLine = line.strip()
            fields = mainLine.split("\t")
            #Convert from CHR:START-STOP:STRAND to CHR START STOP STRAND
            altInfo = fields[9].replace(':','\t')
            altInfo = altInfo.replace('-','\t')
            outputHandle.write(altInfo)
            outputHandle.write("\n")


def fastaExtraction(subInDelbed, fastaSeq, reference):
    print ("Extraction Fasta")
    '''Use the bedtools to extract a sub fasta from a reference and a BED'''
    extractCommand = "bedtools getfasta -fi " + reference + " -bed  " + subInDelbed +" -name+ >" + fastaSeq
    launchingShell(extractCommand)

def BlastOnReference(fastaInput,TEdb,outputBlast):
	print ("Blasting...")
	'''Simple BLAST command launcher'''
	blastCommand="blastn -db " + TEdb + " -query " + fastaInput + " -outfmt 6 -out " + outputBlast
	launchingShell(blastCommand)

def TEidentificationFromBlast(blastresults, seqDict, outputFile, minSize):
	print ("Checking if TE are moving")
	'''Use a simple BLAST parsing to identify complete TE new insertion/deletion'''
	outputHandle = open(outputFile, "w")
	outputHandle.write("#TE\tLocation\tPercId\tFragSize\tRefSize\tPercTotal\n")
	with open(blastresults, "r") as inputHandle:
		for line in inputHandle:
			line = line.strip()
			fields = line.split("\t")
			if int(fields[3]) < ((minSize/100) * len(seqDict[fields[1]])):
				#print ("too short for " + fields[1])
				#exit()
				continue
			else:
				outLine = fields[1] + "\t" + fields[0] + "\t" + fields[2] + "\t" + fields[3] + "\t" + str(len(seqDict[fields[1]])) + "\t" + str((int(fields[3])/len(seqDict[fields[1]]))*100) + "\n"
				outputHandle.write(outLine)


##################################################
## Main code
##################################################
if __name__ == "__main__":
	# Parameters recovery
	parser = argparse.ArgumentParser(prog='svTEidentification.py', description='''This Program annotates the putative TE insertions and deletions in an Assemblytics output''')
	parser.add_argument('-v', '--version', action='version', version='You are using %(prog)s version: ' + version, help=\
						'display svTEidentification.py version number and exit')
	filesreq = parser.add_argument_group('Input mandatory infos for running')
	filesreq.add_argument('-i', '--input', metavar="<filename>", required=True, dest = 'inputFile', help = 'Assemblytics BED file')
	filesreq.add_argument('-r', '--reference', metavar="<filename>", required=True, dest = 'reference', help = 'Reference Sequence file, in fasta')
	filesreq.add_argument('-a', '--alternate', metavar="<filename>", required=True, dest = 'alternate', help = 'Alternate sequence, in fasta')
	filesreq.add_argument('-d', '--database', metavar="<filename>", required=True, dest = 'database', help = 'Multifasta file containing the TE sequences')
	filesreq.add_argument('-o', '--out', metavar="<filename>", required=True, dest = 'outputPrefix', help = 'Prefix for the output files')
	parser.add_argument('-s', '--size', metavar="<Percentage>", required=False, default=90, type=int, dest = 'minimalSize', help = 'Minimal percentage of TE size to be analyzed')




	# Check parameters
	args = parser.parse_args()

	#Welcome message
	print("##########################################################")
	print("# 	      svTEidentification (Version " + version + ")	         #")
	print("##########################################################")

	#Window size for scanning
	minimalSize=args.minimalSize
	print ("Minimal percentage of similarity and size is of " + minimalSize + "%")

	# From relative to absolute paths
	inputFile = relativeToAbsolutePath(args.inputFile)
	reference=relativeToAbsolutePath(args.reference)
	alternate=relativeToAbsolutePath(args.alternate)
	outputPrefix = relativeToAbsolutePath(args.outputPrefix)
	database = relativeToAbsolutePath(args.database)

     #Charging info for TE size
	seqDict = SeqIO.to_dict(SeqIO.parse(database, "fasta"))
	#print(len(seqDict["ZAM"]))

    #Creating the sub BED with indel only
	indelExtraction(inputFile, outputPrefix)

    #Deletion step
	print("Working on Deletions...")
	subDelBed = outputPrefix + "_deletion.bed"
	deletionFastaSeq = outputPrefix + "_deletion.fasta"
	fastaExtraction(subDelBed, deletionFastaSeq, reference)
	#Launching blast command on the database
	print("BLAST control on TE reference database for Deletions")
	delBlastOutput = outputPrefix + "_deletions_vsTEdb.bln"
	BlastOnReference(deletionFastaSeq,database,delBlastOutput)
	deletionOutput = outputPrefix + "_deleted_TE.csv"
	TEidentificationFromBlast(delBlastOutput, seqDict, deletionOutput, minimalSize)


	#Insertion step
	print("Working on Insertions...")
	subInsBed = outputPrefix + "_insertion.bed"
	subInsBedAlt = outputPrefix + "_insertion_ALT.bed"
	BEDconverter(subInsBed,subInsBedAlt)
	insertionFastaSeq = outputPrefix + "_insertion.fasta"
	fastaExtraction(subInsBedAlt, insertionFastaSeq, alternate)
	#Launching blast command on the database
	print("BLAST control on TE reference database for Insertions")
	insBlastOutput = outputPrefix + "_insertions_vsTEdb.bln"
	BlastOnReference(insertionFastaSeq, database, insBlastOutput)
	insertionOutput = outputPrefix + "_insertion_TE.csv"
	TEidentificationFromBlast(insBlastOutput, seqDict, insertionOutput, minimalSize)


print("##########################################################")
print("# 	      FINISHED! 	         #")
print("# 	      Please cite xxx 	         #")
print("##########################################################")
