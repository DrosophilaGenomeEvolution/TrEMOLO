#!/usr/bin/env python3
# -*- coding: utf-8 -*-

######################################################################################################################
#
# Copyright 2019-2020 IRD-CNRS-Lyon1 University
#
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program; if not, see <http://www.gnu.org/licenses/> or
# write to the Free Software Foundation, Inc.,
# 51 Franklin Street, Fifth Floor, Boston,
# MA 02110-1301, USA.
#
# You should have received a copy of the CeCILL-C license with this program.
# If not see <http://www.cecill.info/licences/Licence_CeCILL-C_V1-en.txt>
#
# Intellectual property belongs to authors and IRD, CNRS, and Lyon 1 University  for all versions
# Version 0.1 written by Francois Sabot
#
######################################################################################################################
"""
    Parsing a mummer output to extract variations
    ==========================
    :author:  François Sabot
    :contact: francois.sabot@ird.fr
    :date: 17/11/2020
    :version: 0.1
    :licence: GPLv3 - CeCILL-C
    Script description
    ------------------
    extractVarFromMummer.py will take a double mummer input (delta and show-coords files) to recreate an
    output BED file containing variation.
    It is widely inspired from the excellent job done by Maria Nattestad in her Assemblytics tools:
    https://github.com/MariaNattestad/Assemblytics
    -------
    >>> extractVarFromMummer.py -d deltaFile -s show-coordsFile -o outputName [-m minimalSize -l maximalSize -i True ]

    Help Programm
    -------------
    information arguments:
        - \-h, --help
                        show this help message and exit
        - \-v, --version
                        display svTEidentification.py version number and exit
    Input mandatory infos for running:
        - \- d <filename>, --delta <filename>
                        The delta file issued from Mummer v4.0
        - \- s <filename>, --showcoords <filename>
                        The show-coords file issued from Mummer v4.0 (correspond to the delta file)
        - \-o <filename>, --out <filename>
                        Prefix of output files
    Optional informations
        - \-m <minimalSize>, --minimum <minimalSize>
                        Minimal size for a variation to be recovered  in base pairs (default: 100)
        - \-l <maximalSize>, --maximal <maximalSize>
                        Maximal size for a variation to be recovered  in base pairs - 0 for no limit (default: 0)
        - \-i <True/False>, --interchromosomal <True/False>
                        Showing the interchromosomal variations, True or False (default: True)

"""

# Import
import os

current_dir = os.path.dirname(os.path.abspath(__file__)) + "/"

import subprocess as sp

## Python modules
import argparse
import logging
import gzip


##################################################
# Global variables
##################################################
version = "0.1"
VERSION_DATE = '17/11/2020'
debug = "False"
toolName = "extractVarFromMummer.py"
# debug = "True"

##################################################
# Functions
##################################################

def relativeToAbsolutePath(relative):
    if relative[0] != "/":
        # The relative path is a relative path, ie do not starts with /
        command = "readlink -m " + relative
        absolutepath = sp.check_output(command, shell=True).decode("utf-8").rstrip()
        return absolutepath
    else:
        # Relative is in fact an absolute path, send a warning
        # Relative is in fact an absolute path, send a warning
        absolutepath = relative
        return absolutepath

def str2bool(v):
    if isinstance(v, bool):
       return v
    if v.lower() in ('yes', 'true', 't', 'y', '1'):
        return True
    elif v.lower() in ('no', 'none', 'false', 'f', 'n', '0'):
        return False
    else:
        raise argparse.ArgumentTypeError('Boolean value expected.')

def launchingshell(command):
    # Function to launch and display command
    """Launcher for command"""
    process = sp.Popen(command, stdout=sp.PIPE, shell=True, stderr=sp.PIPE)
    (outputCommand, errorCommand) = process.communicate()
    return outputCommand, errorCommand

def uniqAnchor(delta,minimum,maximum,inter,prefix):
    """Lanch the assemblytics uniq anchor script modified"""
    outputtab = prefix + "/uniq"
    commanduniq="scripts/Assemblytics_uniq_anchor.py --keep-small-uniques --delta " + delta \
    + " --out " + outputtab + " --unique-length " + minimum
    (standardout, errorout) = launchingshell(commanduniq)
    logging.info(standardout)
    if errorout:
        logging.error(errorout)

def betweenAlignment(coords, minimum, maximum, inter, prefix):
    """Take a tabulated output file and create a temporary bed file withe the variation
    Based on the perl script from Assemblytics"""

    tempbed = prefix + "/between.bed"
    outputhandle = open(tempBed, "w")

    try:
        inputhandle = gzip.open(coords, "r")
        logging.info("Opening gzipped file " + coords )
    except:
        inputhandle = open(coords, "r")
        logging.info("Opening plain file " + coords)
    # removing headers
    logging.info("Reading...")
    header1 = inputhandle.readline()
    header2 = inputhandle.readline()
    header3 = inputhandle.readline()
    header4 = inputhandle.readline()
    logging.info("Headers for " + coords + " are:\n\t"
                 + header1 + "\t"
                 + header2 + "\t"
                 + header3 + "\t"
                 + header4)
    previousline = inputhandle.readline().strip()
    previousfields = previousline.split()
    farestend = int(previousfields[1])
    for line in inputhandle:
        fields = line.strip().split()
        if farestend > int(fields[0]):
            # overlap/contained in, next
            previousfields = fields
            # recompute the end of overlap
            farestend = max(farestend,int(fields[1]))
            continue

        localdistref = abs(int(fields[1]) - farestend) + 1
        localdistalt = (int(fields[1]) - int(previousfields[0])) + 1

        if localdistref > minimum and localdistref < maximum and maximum != 0:
            if localdistref > 0:
                if localdistalt > 0:
                    type = "FF"
                elif localdistalt < 0:
                    type = "FR"




        previousfields = fields



        # Convert from CHR:START-STOP:STRAND to CHR START STOP STRAND
        altInfo = fields[9].replace(':', '\t')
        altInfo = altInfo.replace('-', '\t')
        outputHandle.write(altInfo)
        outputHandle.write("\n")

def withinAlignment(delta, minimum, maximum, inter, prefix):

def concatenate(prefix):



# #################################################
# Main code
# #################################################
if __name__ == "__main__":

    # Welcome message
    print("##########################################################")
    print("#       " + toolName + " (Version " + version + ")            #")
    print("##########################################################")

    # Parameters recovery
    parser = argparse.ArgumentParser(prog='toolName',
                                     description='''\nextractVarFromMummer.py will take a double mummer input \
                                     (delta and show-coords files) to recreate an \
                                     output BED file containing variation''')
    parser.add_argument('-v', '--version', action='version',
                        version='You are using ' + toolName + ' version: ' + version,
                        help='display ' + toolName + ' version number and exit')
    # Mandatory ones
    filesreq = parser.add_argument_group('Input mandatory information for running')
    filesreq.add_argument('-d', '--delta', metavar="<filename>", required=True, dest='deltaFile',
                          help='Delta file coming from Mummer 4.0+')
    filesreq.add_argument('-s', '--showcoords', metavar="<filename>", required=True, dest='showCoords',
                          help='Show-coords file coming from Mummer 4.0+, corresponding to the delta file')
    filesreq.add_argument('-o', '--out', metavar="<filename>", required=True,
                          dest='outputPrefix', help='Prefix for the output files')

    # Optional ones
    parser.add_argument('-m', '--minimal', metavar="<minimal>", required=False, default=100, type=int,
                        dest='minimalSize',
                        help='Minimal size for a variation to be analyzed '
                             'default 100')
    parser.add_argument('-l', '--maximal', metavar="<maximal>", required=False, default=0, type=int,
                        dest='maximalSize',
                        help='Maximal size for a variation to be analyzed '
                             '0 = no limit, default')
    parser.add_argument('-i', '--interchromosomal', metavar="<minimal>", required=False, default=True, type=str2bool,
                        dest='interchromosomal',
                        help='Interchromosomal variants to be shown'
                             'default True')
    args = parser.parse_args()

    # From relative to absolute paths
    deltaFile = relativeToAbsolutePath(args.deltaFile)
    showCoords = relativeToAbsolutePath(args.showCoords)
    outputPrefix = relativeToAbsolutePath(args.outputPrefix)

    # Logging and output system
    os.makedirs(outputPrefix, exist_ok=True)
    outputLog = outputPrefix + "/extract.log"
    outputFile = outputPrefix + "/extract.bed"
    logging.basicConfig(filename=outputLog, level=logging.DEBUG, format='%(asctime)s %(message)s')
    logging.info(toolName + " version " + version)
    logging.info("Working on " + deltaFile + " and " + showCoords + " files ")
    logging.info("Output files will be in the " + outputPrefix + "folder")


    # Parameters to conserve variants
    minimalSize = args.minimalSize
    maximalSize = args.maximalSize
    if minimalSize > maximalSize > 0:
        print("Minimal size must be higher than maximal one!  \n\n\t\tExiting")
        logging.error("Minimal size must be higher than maximal one!  \n\n\t\tExiting")

    logging.info("Minimal size is of " + str(minimalSize) + "bp")
    if maximalSize == 0:
        logging.info("No limit to maximal size")
    else:
        logging.info("Maximal size is of " + str(maximalSize) + "bp")

    interchromosomal = args.interchromosomal
    if interchromosomal:
        logging.info("Interchromosomal variants will be conserved in the output file")
    else:
        logging.info("Interchromosomal variants will not be conserved in the output file")

    # The following part is mostly inspired from Assemblytics but adapted to the show-coords file and use their scripts
    # converting in a tabular file
    uniqAnchor(deltaFile,minimalSize,maximalSize,interchromosomal, outputPrefix)

    exit(0)
    # Starting between alignment analysis
    betweenAlignment(showCoords, minimalSize, maximalSize, interchromosomal, outputPrefix)
    # Continuing within alignment analysis
    withinAlignment(deltaFile, minimalSize, maximalSize, interchromosomal, outputPrefix)
    # Mixing the two
    concatenate(outputPrefix)
