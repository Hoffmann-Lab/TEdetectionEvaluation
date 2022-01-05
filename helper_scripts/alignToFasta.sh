#!/usr/bin/env bash
alignFile=$1
referenceGenome=$2

echo given align file: $alignFile
echo given reference genome: $referenceGenome

# translate align file to bed and gtf
python3 align_parser.py -a $alignFile


# Build file prefix
filePrefix=${alignFile%%.*}

# bed file to fasta file
# The awk command removes the position in name added by bedtools
echo bedtools getfasta -fi $referenceGenome -bed $filePrefix.repeatAnnotation.bed -name -fo $filePrefix.repeats.fa
bedtools getfasta -fi $referenceGenome -bed $filePrefix.repeatAnnotation.bed -name | \
    awk 'BEGIN{FS="::"}{if($1~">"){printf $1"\n"}else{print $0}}' > $filePrefix.repeats.fa
