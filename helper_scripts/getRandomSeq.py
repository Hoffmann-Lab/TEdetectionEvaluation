#!/gsc/biosw/bin/python3.7

from pylibs.seq_handling import load_seq_file
import sys
import random
import argparse
from Bio import SeqIO


def parse_options():

    parser = argparse.ArgumentParser(prog="reext", description="REpetitive\
            element EXtractor (reext) is a tool to extract the reference\
            sequences of repetitive elements annotated by repeatmasker.")

    parser.add_argument("-f", type=str, dest="fastaFile")
    parser.add_argument("-n", type=int, dest="number")
    parser.add_argument("-l", type=int, default = 100, dest="minLength")
    parser.add_argument("-o", type=str, default = "randomSet", dest="filePrefix")
    parser.add_argument("-s", type=int, default = -1, dest="seed")

    args = parser.parse_args()

    return parser, args


parser, args = parse_options()

def set_seed(seed):

    if seed == -1:

        seed = random.choice(list(range(1,1000)))
        print("Seed: {}".format(seed))
    random.seed(seed)


#set seed
set_seed(args.seed)

fastaEntries = load_seq_file(args.fastaFile,args.minLength)

if fastaEntries == None:
    sys.exit("No sequence file")

if len(fastaEntries) < args.number:
    sys.exit("There are less entries in your reference file as the requested \
sample number.")

randomSet = random.sample(fastaEntries, args.number)

assert len(randomSet) == args.number, "ooops the random set is smaller than \
the requested number"

with open("{}.{}.mi{}bp.fa".format(args.filePrefix, args.number,
    args.minLength), "w") as output_handle:

    SeqIO.write(randomSet, output_handle, "fasta") 
