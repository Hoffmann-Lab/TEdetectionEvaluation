#!/gsc/biosw/bin/python3.7

########################
#
# author: Robert Schwarz
# mail: rschwarz@leibniz-fli.de
# initial date:  Fri Oct 10 14:40:18 CEST 2019
#
# description:  Reext - repetitive element extraction, is a script that 
#               uses the information from an .align-file (output Repeatmasker)
#               to extract the sequences of the repetitive elements and their
#               kimura distance. The .align file is transcribed into
#               a .bed -file, which is used by bedtools to extract
#               the sequences of each TE instance from the reference genome
#               and stores them in a new fasta file. The name of the instances
#               is a combination of a lot of information like, location, family,
#               kimura distance to name some examples.
# 
#
#
# .align File of Repeatmasker
#   This file contains, in contrast to the .out file, the alignment and
#   the kimura distance.
#
#
#########
# To-Do #
#########
#
# - if a bed file for all TE is existing in the output directory avoid to
#   create a new one
# - add the possibility to read the .out-file
# - add the possibility to determine the percentage that should be
#   simulated
# - store the information if a TE instances are not found in the 
#   reference genome. This information is called bei bedtools. This
#   leads to the fact that a few instances less are stored in 
#   the reference file
# - write a method that checks if the arguments are declared that are needed
# - fill the args and argv into the info dict so that the log method neads only one
#   attribute

import sys
import os
import argparse

import pylibs.fileparser as fparse
#import pylibs.administration as admin



def parse_options():

    parser = argparse.ArgumentParser(prog="reext", description="REpetitive\
            element EXtractor (reext) is a tool to extract the reference\
            sequences of repetitive elements annotated by repeatmasker.")

    parser.add_argument("-a", type=str, dest="align_file")

    parser.add_argument("-o", type=str, default="result", dest="output")

    parser.add_argument("-mi", type=int, default = 0, dest="min_length")

    args = parser.parse_args()

    return parser, args


def translate_to_annotation_files(alignEntries, filePrefix):

   
    bedFileName = "{}.repeatAnnotation.bed".format(filePrefix)
    gtfFileName = "{}.repeatAnnotation.gtf".format(filePrefix)
   
    with open(bedFileName, "w") as bed_handle, \
    open(gtfFileName, "w") as gtf_handle:
        
        for instance in alignEntries:
            
            bedLine = alignEntries[instance].translate_to_bed()
            gtfLine = alignEntries[instance].translate_to_gtf()
            
            gtf_handle.write(gtfLine + "\n")
            bed_handle.write(bedLine + "\n")

def load_align(align_file):
        
        align_entries = {}
        
        print("load align file...")

        for instance in fparse.parse_align_file(align_file):
            
            if instance.seq_id in align_entries:

                align_entries[instance.seq_id].occurence += 1

            else:
                align_entries[instance.seq_id] = instance
        
        print("Done.")
        return align_entries

def main():

    parser, args = parse_options()
    
    filePrefix = args.align_file.split(".")[0]

    alignEntries = load_align(args.align_file)
    translate_to_annotation_files(alignEntries, filePrefix)

if __name__ == "__main__":
    main()
