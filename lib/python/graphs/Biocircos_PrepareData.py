#!/usr/bin/python 
""" Module/Script Description

This Module/Script can help users prepare data for biocircos.

Copyright (c) 2015
This code is free software; you can redistribute it and/or modify it
under the terms of the BSD License (see the file COPYING included with
the distribution).
@status:  experimental
@version: $Revision$
"""

# ------------------------------------
# python modules
# ------------------------------------

import sys
import os

# ------------------------------------
# constants
# ------------------------------------

# ------------------------------------
# Misc functions
# ------------------------------------

# ------------------------------------
# Classes
# ------------------------------------

# ------------------------------------
# Main
# ------------------------------------

if __name__=="__main__":
    if len(sys.argv)==1:
        print( "Usage: "+sys.argv[0]+" SCATTER/SNP/ARC/HEATMAP/LINE/HISTOGRAM/LINK/CNV  SCATTER01.txt/SNP01.txt/ARC01.txt/HEATMAP01.txt/LINE01.txt/HISTOGRAM01.txt/LINK01.txt/CNV01.txt > SCATTER01.js/SNP01.js/ARC01.js/HEATMAP01.js/LINE01.js/HISTOGRAM01.js/LINK01.js/CNV01.js")
    else:
        datatype=sys.argv[1]
        dataname=sys.argv[2]

        # Prepare data for SCATTER demo (gallery04_SCATTER01.js).
        if datatype=="SCATTER":
            print( "var "+dataname.split(".")[0]+" = [ \""+dataname.split(".")[0]+"\" , {")
            print( "  SCATTERRadius: 220,")
            print( "  innerCircleSize: 1,")
            print( "  outerCircleSize: 7,")
            print( "  innerCircleColor: \"red\",")
            print( "  outerCircleColor: \"#CC3399\",")
            print( "  innerPointType: \"circle\", //circle,rect")
            print( "  outerPointType: \"circle\", //circle,rect")
            print( "  innerrectWidth: 2,")
            print( "  innerrectHeight: 2,")
            print( "  outerrectWidth: 10,")
            print( "  outerrectHeight: 10,")
            print( "  outerCircleOpacity: 1,")
            print( "  random_data: 0,")
            print( "  SCATTERAnimationDisplay: false,")
            print( "  SCATTERAnimationTime: 2000,")
            print( "  SCATTERAnimationDelay: 20,")
            print( "  SCATTERAnimationType: \"bounce\"")
            print( "} , [")

            fh=open(sys.argv[2])
            for line in fh:
                line=line.rstrip(os.linesep)
                if line[0]!="#": #Title line started with "#"
                    line=line.split("\t")
                    print( "  {chr: \""+line[0]+"\", start: \""+line[1]+"\", end: \""+line[2]+"\", name: \""+line[3]+"\", des: \""+line[4]+"\"},")
            print( "]];")

        # Prepare data for SNP demo (gallery05_SNP01.js).
        if datatype=="SNP":
            print( "var "+dataname.split(".")[0]+" = [ \""+dataname.split(".")[0]+"\" , {")
            print( "  maxRadius: 205,")
            print( "  minRadius: 173,")
            print( "  SNPFillColor: \"#9400D3\",")
            print( "  PointType: \"circle\",")
            print( "  circleSize: 2,")
            print( "  rectWidth: 2,")
            print( "  rectHeight: 2,")
            print( "  SNPAnimationDisplay: false,")
            print( "  SNPAnimationTime: 2000,")
            print( "  SNPAnimationDelay: 20,")
            print( "  SNPAnimationType: \"bounce\"")
            print( "} , [")

            fh=open(sys.argv[2])
            colorDict={0:"rgb(153,102,0)", 1:"rgb(102,102,0)", 2:"rgb(153,153,30)", 3:"rgb(204,0,0)", 4:"rgb(255,0,0)", 5:"rgb(255,0,204)", 6:"rgb(255,204,204)", 7:"rgb(255,153,0)", 8:"rgb(255,204,0)", 9:"rgb(255,255,0)", 10:"rgb(204,255,0)", 11:"rgb(0,255,0)", 12:"rgb(53,128,0)", 13:"rgb(0,0,204)", 14:"rgb(102,153,255)", 15:"rgb(153,204,255)", 16:"rgb(0,255,255)", 17:"rgb(204,255,255)", 18:"rgb(153,0,204)", 19:"rgb(204,51,255)", 20:"rgb(204,153,255)", 21:"rgb(102,102,102)", 22:"rgb(153,153,153)", 23:"rgb(204,204,204)"}
            chrList=[]
            for line in fh:
                line=line.rstrip(os.linesep)
                if line[0]!="#":#Title line started with "#"
                    line=line.split("\t")
                    if(line[0] in chrList):
                        pass
                    else:
                        chrList.append(line[0])
                    print( "  {chr: \""+line[0]+"\", pos: \""+line[1]+"\", value: \""+line[2]+"\", des: \""+line[3]+"\", color: \""+colorDict[chrList.index(line[0])]+"\"},")
            print( "]];")

        # Prepare data for ARC demo (gallery06_ARC01.js).
        if datatype=="ARC":
            print( "var "+dataname.split(".")[0]+" = [ \""+dataname.split(".")[0]+"\" , {")
            print( "  innerRadius: -55,")
            print( "  outerRadius: -45,")
            print( "} , [")

            fh=open(sys.argv[2])
            for line in fh:
                line=line.rstrip(os.linesep)
                if line[0]!="#":#Title line started with "#"
                    line=line.split("\t")
                    print( "  {chr: \""+line[0]+"\", start: \""+line[1]+"\", end: \""+line[2]+"\", color: \""+line[3]+"\", des: \""+line[4]+"\"},")
            print( "]];")

        # Prepare data for HEATMAP demo (gallery07_HEATMAP01.js).
        if datatype=="HEATMAP":
            print( "var "+dataname.split(".")[0]+" = [ \""+dataname.split(".")[0]+"\" , {")
            print( "  innerRadius: -25,")
            print( "  outerRadius: -65,")
            print( "  maxColor: \"red\",")
            print( "  minColor: \"yellow\"")
            print( "} , [")

            fh=open(sys.argv[2])
            for line in fh:
                line=line.rstrip(os.linesep)
                if line[0]!="#": #Title line started with "#"
                    line=line.split("\t")
                    print( "  {chr: \""+line[0]+"\", start: \""+line[1]+"\", end: \""+line[2]+"\", name: \""+line[3]+"\", value: \""+line[4]+"\"},")
            print( "]];")

        # Prepare data for LINE demo (gallery08_LINE01.js).
        if datatype=="LINE":
            print( "var "+dataname.split(".")[0]+" = [ \""+dataname.split(".")[0]+"\" , {")
            print( "  maxRadius: 220,")
            print( "  minRadius: 170,")
            print( "  LineColor: \"#EEAD0E\",")
            print( "  LineWidth: 2,")
            print( "} , [")

            fh=open(sys.argv[2])
            for line in fh:
                line=line.rstrip(os.linesep)
                if line[0]!="#": #Title line started with "#"
                    line=line.split("\t")
                    print( "  {chr: \""+line[0]+"\", pos: \""+line[1]+"\", value: \""+line[2]+"\"},")
            print( "]];")
        
        # Prepare data for HISTOGRAM demo (gallery09_HISTOGRAM.js).
        if datatype=="HISTOGRAM":
            print( "var "+dataname.split(".")[0]+" = [ \""+dataname.split(".")[0]+"\" , {")
            print( "  maxRadius: 240,")
            print( "  minRadius: 205,")
            print( "  histogramFillColor: \"#FF6666\",")
            print( "} , [")

            fh=open(sys.argv[2])
            for line in fh:
                line=line.rstrip(os.linesep)
                if line[0]!="#": #Title line started with "#"
                    line=line.split("\t")
                    print( "  {chr: \""+line[0]+"\", start: \""+line[1]+"\", end: \""+line[2]+"\", name: \""+line[3]+"\", value: \""+line[4]+"\"},")
            print( "]];")
        
        # Prepare data for LINK demo (gallery10_LINK01.js).
        if datatype=="LINK":
            print( "var "+dataname.split(".")[0]+" = [ \""+dataname.split(".")[0]+"\" , {")
            print( "  LinkRadius: 140,")
            print( "  LinkFillColor: \"#F26223\",")
            print( "  LinkWidth: 3,")
            print( "  displayLinkAxis: true,")
            print( "  LinkAxisColor: \"#B8B8B8\",")
            print( "  LinkAxisWidth: 0.5,")
            print( "  LinkAxisPad: 3,")
            print( "  displayLinkLabel: true,")
            print( "  LinkLabelColor: \"red\",")
            print( "  LinkLabelSize: 13,")
            print( "  LinkLabelPad: 8,")
            print( "} , [")

            fh=open(sys.argv[2])
            for line in fh:
                line=line.rstrip(os.linesep)
                if line[0]!="#": #Title line started with "#"
                    line=line.split("\t")
                    print( "  {fusion: \""+line[0]+"\", g1chr: \""+line[1]+"\", g1start: \""+line[2]+"\", g1end: \""+line[3]+"\", g1name: \""+line[4]+"\", g2chr: \""+line[5]+"\", g2start: \""+line[6]+"\", g2end: \""+line[7]+"\", g2name: \""+line[8]+"\"},")
            print( "]];")

        # Prepare data for CNV demo (Figure_CNV_CNV.js).
        if datatype=="CNV":
            print( "var "+dataname.split(".")[0]+" = [ \""+dataname.split(".")[0]+"\" , {")
            print( "  maxRadius: 155,")
            print( "  minRadius: 116,")
            print( "  CNVwidth: 2,")
            print( "  CNVColor: \"#4876FF\",")
            print( "} , [")

            fh=open(sys.argv[2])
            for line in fh:
                line=line.rstrip(os.linesep)
                if line[0]!="#": #Title line started with "#"
                    line=line.split("\t")
                    print( "  {chr: \""+line[0]+"\", start: \""+line[1]+"\", end: \""+line[2]+"\", value: \""+line[3]+"\"},")
            print( "]];")
