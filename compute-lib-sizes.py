#!/usr/bin/env python
#=========================================================================
# This is OPEN SOURCE SOFTWARE governed by the Gnu General Public
# License (GPL) version 3, as described at www.opensource.org.
# Copyright (C)2021 William H. Majoros <bmajoros@alumni.duke.edu>
#=========================================================================
from __future__ import (absolute_import, division, print_function, 
   unicode_literals, generators, nested_scopes, with_statement)
from builtins import (bytes, dict, int, list, object, range, str, ascii,
   chr, hex, input, next, oct, open, pow, round, super, filter, map, zip)
# The above imports should allow this program to run in both Python 2 and
# Python 3.  You might need to update your version of module "future".
import sys
import ProgramName
import gzip
from Rex import Rex
rex=Rex()

#=========================================================================
# main()
#=========================================================================
if(len(sys.argv)!=2):
    exit(ProgramName.get()+" <in-counts.txt.gz>\n")
(infile,)=sys.argv[1:]

dnaSums=None; rnaSums=None
with gzip.open(infile,"rt") as IN:
    header=IN.readline()
    if(not rex.find("DNA=(\d+)\s+RNA=(\d+)",header)):
        raise Exception("expecting RNA=# DNA=# header")
    numDna=int(rex[1]); numRna=int(rex[2])
    dnaSums=numDna*[0]; rnaSums=numRna*[0]
    for line in IN:
        fields=line.rstrip().split()
        dnaCounts=[int(x) for x in fields[:numDna]]
        rnaCounts=[int(x) for x in fields[numDna:(numDna+numRna)]]
        for i in range(numDna): dnaSums[i]+=dnaCounts[i]
        for i in range(numRna): rnaSums[i]+=rnaCounts[i]
for i in range(numDna):
    label="input_rep"+str(i+1)
    print(label,dnaSums[i],sep="\t")
for i in range(numRna):
    label="output_rep"+str(i+1)
    print(label,rnaSums[i],sep="\t")


