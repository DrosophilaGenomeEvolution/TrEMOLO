#!/usr/bin python3
# -*- coding: utf-8 -*-

###################################################################################################################################
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
#If not see <http://www.cecill.info/licences/Licence_CeCILL-C_V1-en.txt>
#
# Intellectual property belongs to authors and IRD, CNRS, and Lyon 1 University  for all versions
# Version 0.1 written by Mourdas Mohamed
#                                                                                                                                   
####################################################################################################################################

"""
    CONVERT FASTQ TO FASTA
    ==========================
    :author:  Mourdas MOHAMED 
    :contact: mourdas.mohamed@igh.cnrs.fr
    :date: 01/06/2020
    :version: 0.1
    Script description
    ------------------
    fastq_to_fasta.py convert fastq to fasta
    -------
    >>> fastq_to_fasta.py file.fastq file_out.fasta
    
    Help Programm
    -------------
    usage: fastq_to_fasta.py [-h] fastq_file fasta_out_file

    Convert fastq to fasta

    positional arguments:
      fastq_file      input fastq file
      fasta_out_file  fasta output file

    optional arguments:
      -h, --help      show this help message and exit
"""


import sys
import os
import pandas
import re
import argparse

#parse args
parser = argparse.ArgumentParser(description="Convert fastq to fasta")

#MAIN ARGS
parser.add_argument("fastq_file", type=str,
                    help="input fastq file")
parser.add_argument("fasta_out_file", type=str,
                    help="fasta output file")

args = parser.parse_args()

name_fastq = args.fastq_file
name_fasta = args.fasta_out_file
file_out   = open(name_fasta, "w")
file       = open(name_fastq, "r")
line       = file.readline()

sequence          = ""
qual              = ""
while line:
    if line[0] == "@":
        sequence = ""
        qual     = ""
        file_out.write(">"+line[1:])
        line = file.readline()
        while line and line[0] != "+":
            sequence += line.replace("\n", "")
            line = file.readline()
    
    if line and line[0] == "+":
        file_out.write(sequence+"\n")
        line = file.readline()
        while line and len(qual) < len(sequence):
            qual += line.replace("\n", "")
            line = file.readline()

file_out.close()
file.close()