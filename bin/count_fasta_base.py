#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
From: https://techoverflow.net/2013/10/24/a-simple-tool-for-fasta-statistics/
By: Ulrich Josef Stefan Köhler
fasta-stats.py: Utility script to count number of nucleotides/Aminoacids in a FASTA files.
Changelog:
    1.1: Python3-ready
"""
from __future__ import with_statement
import sys
import argparse
import gzip

#Counter is used for per-character statistics
from collections import Counter
__author__    = "Uli Koehler & Anton Smirnov"
__copyright__ = "Copyright 2013 Uli Koehler & Anton Smirnov"
__license__   = "Apache v2.0"
__version__   = "1.1"
def printSequenceStats(fileName, sequenceName, charOccurrenceMap, totalCharCount):
    """
    Print details about a sequence to stdout.
    Called from inside parseFile().
    Keyword arguments:
        sequenceName: The name of the sequence to print
        charOccurrenceMap: A dict-like object that contains a char -->
                            number of occurrences mapping
        totalCharCount: The overall character count of the sequence
    """
    # print ("There are {} characters in total:".format(totalCharCount))
    
    relativeFrequencies = []
    for char in sorted(charOccurrenceMap.keys()):
        charCount = charOccurrenceMap[char]
        relativeFrequency = charCount * 1.0 / totalCharCount
        relativeFrequencies.append(relativeFrequency)
        print ("{}".format(relativeFrequency))
    #For nucleotide sequences (ATGC only), also print A+T vs G+C count
    # if sorted(charOccurrenceMap.keys()) == ["A","C","G","T"]:
    #     #Print A+T count
    #     atCount = charOccurrenceMap["A"] + charOccurrenceMap["T"]
    #     atRelFrequency = atCount * 100.0 / totalCharCount
    #     print ("\tA+T : {} = {}%".format(atCount, atRelFrequency))
    #     #Print G+C count
    #     gcCount = charOccurrenceMap["G"] + charOccurrenceMap["C"]
    #     atRelFrequency = gcCount * 100.0 / totalCharCount
    #     print ("\tG+C : {} = {}%%".format(gcCount, atRelFrequency))
def parseFile(filename, caseSensitive=False, charWhitelist=None):
    """
    Parse a FASTA fil and call printRe
    """
    #Set to the header line, with ">" removed from the beginning
    sequenceName = None
    #Key: character, value: Number of occurrences
    charOccurrenceMap = Counter()
    #The number of characters in the current sequence, not including \n
    charCount = 0
    #Keep track of consecutive comments, because they are appended
    previousLineWasComment = False
    #Open and iterate the file, auto- detect gzip
    openFunc = gzip.open if filename.endswith(".gz") else open
    with openFunc(filename, "r") as infile:
        for line in infile:
            line = line.strip()
            #Be super-compatible with the original specification
            if line.startswith(">") or line.startswith(";"):
                #Process previous sequence, if any
                if sequenceName is not None:
                    printSequenceStats(filename, sequenceName, charOccurrenceMap, charCount)
                    charOccurrenceMap = Counter()
                    charCount = 0
                #Take the entire comment line as (new) sequence ID (with ">" stripped)
                #Concatenate consecutive sequence lines
                if previousLineWasComment: #Append -- add one space between to normalize whitespace count
                    sequenceName += " " + line[1:].strip()
                else:
                    sequenceName = line[1:].strip()
                previousLineWasComment = True
            else: #Line belongs to the sequence
                previousLineWasComment = False
                #Line has been stripped before, so we can count directly
                #Increment per-character stats (character occurrences)
                for char in line:
                    #Skip any character not in the whitelist, if whitelist (--only) is enabled
                    if charWhitelist is not None and not char in charWhitelist:
                        continue
                    #We can only count after whitelite filter
                    charCount += 1
                    #In case-insensitive mode (default) count uppercased chars only
                    char = char if caseSensitive else char.upper()
                    charOccurrenceMap[char] += 1
        #The last line has been read, print the last record, if any
        if sequenceName is not None:
            printSequenceStats(filename, sequenceName, charOccurrenceMap, charCount)
if __name__ == "__main__":
    #Allow single or multiple files to be specified
    parser = argparse.ArgumentParser(description='Compute simple statistics for FASTA files.')
    parser.add_argument('infiles', nargs='+', help='FASTA files (.fa, .fa.gz) to generate statistics for')
    parser.add_argument('--case-sensitive', action='store_true', help='Count characters in a case-sensitive way. Disabled per default.')
    parser.add_argument('-o','--only', help='If this option is supplied (e.g. set to \'ATGC\'), characters not in the set will be ignored for all statistics')
    args = parser.parse_args()
    #Process all FASTA files
    for infile in args.infiles:
        parseFile(infile, caseSensitive=args.case_sensitive, charWhitelist=args.only)
