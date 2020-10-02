#!/usr/bin/python2.7

# Script by Devon Ryan ( https://www.biostars.org/p/95929/ )
# for retaining filtered alignments that consist of both reads in a pair
# It achieves this by checking for matching consecutive read names corresponding to a read pair
  
import csv
import sys

f = csv.reader(sys.stdin, dialect="excel-tab")
of = csv.writer(sys.stdout, dialect="excel-tab")
last_read = None
for line in f :
    # Take care of the header
    if(line[0][0] == "@") :
        of.writerow(line)
        continue

    if(last_read == None) :
        last_read = line
    else :
        if(last_read[0] == line[0]) :
            of.writerow(last_read)
            of.writerow(line)
            last_read = None
        else :
            last_read = line
