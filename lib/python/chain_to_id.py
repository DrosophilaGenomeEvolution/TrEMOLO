#!/home/mourdas/anaconda3/bin/python3
import sys
import os

chain=str(sys.stdin.readlines())

s=0
for index, i in enumerate(chain):
    for e in range(4):
        s += (e+ord(i)+(7+index/(5-ord(i))))

s = ( (s - int(s)) * 100000 ) + int(s)
print(int(s))