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


import os
import sys
import argparse

#our
from utils import *

parser = argparse.ArgumentParser()

parser.add_argument("TSD", type=str,
                    help="fasta file of flank")
parser.add_argument("empty_site", type=str,
                    help="file of sequence TE")
parser.add_argument("pos_svi_init", type=int,
                    help="size of flankan")
parser.add_argument("position_on_genome", type=int,
                    help="id of variant rare")

args  = parser.parse_args()


TSD_POSITION = find_svi_position(args.TSD, args.empty_site, args.pos_svi_init, args.position_on_genome)
print("TSD_POSITION:" + str(TSD_POSITION) + "\n\n")