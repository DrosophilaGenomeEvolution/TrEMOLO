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
    >>> extractVarFromMummer.py -d deltaFile -o outputName [-m minimalSize -l maximalSize -i True ]

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
import subprocess as sp
import sys
import argparse
import logging
import gzip
#import yoda_powers

current_dir = os.path.dirname(os.path.abspath(__file__)) + "/"
softwarelocation = pathname = os.path.dirname(sys.argv[0])

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
    commanduniq=softwarelocation + "/scripts/Assemblytics_uniq_anchor.py --keep-small-uniques --delta " + delta + " --out " + outputtab + " --unique-length " + str(minimum)
    # (standardout, errorout) = launchingshell(commanduniq)
    # logging.info(standardout)
    # if errorout:
    #     logging.error(errorout)
    return outputtab + ".coords.tab"
    # TODO be careful, here duplications are not taken into account


def betweenAlignment(coords, minimum, maximum, inter, prefix):
    """Take a tabulated output file and create a temporary bed file withe the variation
    Based on the perl script from Assemblytics"""

    tempbed = prefix + "/between.bed"
    outputhandle = open(tempbed, "w")

    #Gathering infos
    logging.info("Reading...")
    fieldsList=("rstart","rend","qstart","qend","rlen","qlen","rid","qid")
    alignments = dict()
    numalignments = 0
    with open(coords, "r") as input:
        for line in input:
            localdict = dict(zip(fieldsList, line.strip().split()))
            localdict["str"] = line.strip()
            localdict["qidx"] = 0
            localdict["qrc"] = False if localdict["qend"] > localdict["qstart"] else True
            (qid, rid) = (localdict["qid"], localdict["rid"])
            if qid in alignments:
                if rid in alignments[qid]:
                    alignments[qid][rid].append(dict(localdict))
                else:
                    alignments[qid][rid] = [dict(localdict)]
            else:
                alignments[qid] = dict()
                alignments[qid][rid] = [dict(localdict)]
            numalignments += 1
    logging.info("Loaded " + str(numalignments) + " alignments...")

    # Parsing infos
    (candidatefusion, candidatesv, svIdCounter) = (0, 0, 0)
    svstats = dict()
    for qid in sorted(alignments.keys()):
        refs = sorted(alignments[qid].keys())
        qalign = list()
        for rid in refs:
            for aln in alignments[qid][rid]:
                qalign.append(aln)
        qalign = sorted(qalign, key=lambda k: k["qstart"])
        i = 0
        for aln in qalign:
            qalign[i]["qidx"] = i
            i += 1

        # scan for SV
        if len(qalign) > 1:
            # Query has more than one aln, so has SV
            j = 1
            while j < len(qalign):
                previousaln = qalign[j-1]
                currentaln = qalign[j]

                previousline = previousaln["str"]
                currentline = currentaln["str"]

                rid = previousaln["rid"]

                if int(previousaln["rlen"]) >= minimum and int(currentaln["rlen"]) >= minimum:
                    # aln before and after on the reference are longer than the minimal requested here
                    svIdCounter += 1
                    (rpos, qpos, previousposition, currentposition, previousstrand, currentstrand)= ("", "", "", "", "", "")
                    qdist = 0
                    rdist = 0

                    print("+" + previousaln["rlen"] + " " + currentaln["rlen"])
                    print(previousaln["qrc"])
                    print(currentaln["qrc"])
                    if previousaln["qrc"] is False and currentaln["qrc"] is False:
                        svtype = "FF"
                        qdist = int(currentaln["qstart"]) - int(previousaln["qend"])
                        rdist = int(currentaln["rstart"]) - int(previousaln["rend"])
                        print('.')
                        if rdist >= 0:
                            rpos = "%s:%s-%s:+" % (rid, previousaln["rend"], currentaln["rstart"])
                        else:
                            rpos = "%s:%s-%s:-" % (rid, previousaln["rstart"], currentaln["rend"])
                        if qdist >= 0:
                            qpos = "%s:%s-%s:+" % (qid, previousaln["qend"], currentaln["qstart"])
                        else:
                            qpos = "%s:%s-%s:-" % (qid, currentaln["qstart"], previousaln["qend"])
                        previousposition = previousaln["rend"]
                        currentposition = currentaln["rstart"]
                        previousstrand = "+"
                        currentstrand = "-"
                    elif previousaln["qrc"] and currentaln["qrc"]:
                        svtype = "RR"
                        qdist = int(currentaln["qend"]) - int(previousaln["qstart"])
                        rdist = int(previousaln["rstart"]) - int(currentaln["rend"])
                        print('.')
                        if rdist >= 0:
                            rpos = "%s:%s-%s:+" % (rid, previousaln["rend"], currentaln["rstart"])
                        else:
                            rpos = "%s:%s-%s:-" % (rid, previousaln["rstart"], currentaln["rend"])
                        if qdist >= 0:
                            qpos = "%s:%s-%s:+" % (qid, previousaln["qstart"], currentaln["qend"])
                        else:
                            qpos = "%s:%s-%s:-" % (qid, currentaln["qend"], previousaln["qstart"])
                        previousposition = previousaln["rstart"]
                        currentposition = currentaln["rend"]
                        previousstrand = "-"
                        currentstrand = "+"
                    elif previousaln["qrc"] is False and currentaln["qrc"]:
                        svtype = "FR"
                        qdist = int(currentaln["qend"]) - int(previousaln["qend"])
                        rdist = int(currentaln["rstart"]) - int(previousaln["rend"])
                        print('.')
                        if rdist >= 0:
                            rpos = "%s:%s-%s:+" % (rid, previousaln["rend"], currentaln["rstart"])
                        else:
                            rpos = "%s:%s-%s:-" % (rid, currentaln["rstart"], previousaln["rend"])
                        if qdist >= 0:
                            qpos = "%s:%s-%s:+" % (qid, previousaln["qend"], currentaln["qend"])
                        else:
                            qpos = "%s:%s-%s:-" % (qid, currentaln["qend"], previousaln["qstart"])
                        previousposition = previousaln["rend"]
                        currentposition = currentaln["rend"]
                        previousstrand = "+"
                        currentstrand = "+"
                    elif previousaln["qrc"] and currentaln["qrc"] is False:
                        svtype = "RF"
                        qdist = int(previousaln["qend"]) - int(currentaln["qend"])
                        rdist = int(currentaln["rstart"]) - int(previousaln["rend"])
                        print('.')
                        if rdist >= 0:
                            rpos = "%s:%s-%s:+" % (rid, previousaln["rend"], currentaln["rstart"])
                        else:
                            rpos = "%s:%s-%s:-" % (rid, currentaln["rstart"], previousaln["rend"])
                        if qdist >= 0:
                            qpos = "%s:%s-%s:+" % (qid, currentaln["qend"], previousaln["qend"])
                        else:
                            qpos = "%s:%s-%s:-" % (qid, previousaln["qend"], currentaln["qstart"])
                        previousposition = previousaln["rstart"]
                        currentposition = currentaln["rstart"]
                        previousstrand = "-"
                        currentstrand = "-"
                    else:
                        logging.error("ERROR: unknown SV: \n\t\t" + previousline + "\n\t\t" + currentline)

                    totaldist = rdist + qdist
                    typeguess = ""
                    absEventSize = abs(rdist - qdist)

                    if previousaln["rid"] is not currentaln["rid"]:
                        typeguess = "Interchromosomal"
                        rdist = 0
                    else:
                        if previousaln["str"] is currentaln["str"]:
                            typeguess = "Inversion"
                            absEventSize = rdist
                        elif qdist > rdist:
                            if -1 * minimum < rdist < minimum and qdist > -1 * minimum:
                                typeguess = "Insertion"
                            else:
                                if rdist <0 or qdist <0:
                                    typeguess = "TandemExpansion"
                                else:
                                    typeguess = "RepeatExpansion"
                        elif qdist < rdist:
                            if rdist > -1 * minimum and -1 * minimum < qdist < minimum:
                                typeguess = "Deletion"
                            else:
                                if rdist < 0 or qdist < 0:
                                    typeguess = "TandemContraction"
                                else:
                                    type = "RepeatContraction"
                        else:
                            typeguess = "None"
                        if absEventSize > maximum:
                            typeguess = "Longrange"
                            # TODO RECHECK HERE: what about long range ?

                    if typeguess is not "Inversion" and typeguess is not "Interchromosomal" and typeguess is not "None" and absEventSize >= minimum:
                        refstart = min(int(previousaln["rstart"]), int(currentaln["rstart"]))
                        refend = max(int(previousaln["rend"]), int(currentaln["rend"]))
                        if refstart == refend:
                            refend = refstart + 1
                        outlist = (currentaln["rid"], refstart, refend, "Assemblytics_var_" + str(svIdCounter), absEventSize, typeguess, rdist, qdist, qpos)
                        outputhandle.write("\t".join(outlist) + "\n")
                    elif typeguess is "Inversion":
                        # TODO deals with inversion

                    elif typeguess is "Interchromosomal":
                        # TODO deals with interchromosomal

                    candidatesv += 1
                j += 1



    



def withinAlignment(delta, minimum, maximum, inter, prefix):
    return

def concatenate(prefix):
    return



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
    outputPrefix = relativeToAbsolutePath(args.outputPrefix)

    # Logging and output system
    os.makedirs(outputPrefix, exist_ok=True)
    outputLog = outputPrefix + "/extract.log"
    outputFile = outputPrefix + "/extract.bed"
    logging.basicConfig(filename=outputLog, level=logging.DEBUG, format='%(asctime)s %(message)s')
    logging.info(toolName + " version " + version)
    logging.info("Working on " + deltaFile + " file ")
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
    tabularfile = uniqAnchor(deltaFile,minimalSize,maximalSize,interchromosomal, outputPrefix)

    # Starting between alignment analysis
    betweenAlignment(tabularfile, minimalSize, maximalSize, interchromosomal, outputPrefix)
    # Continuing within alignment analysis
    withinAlignment(deltaFile, minimalSize, maximalSize, interchromosomal, outputPrefix)
    # Mixing the two
    concatenate(outputPrefix)
