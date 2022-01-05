########################
#
# author: Robert Schwarz
# mail: rschwarz@leibniz-fli.de
# initial date:  Fri Oct 8 14:40:18 CEST 2019
#
# description: This is a lib to parse annotation files.
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

# - class out_entry !?
# - align_entry can then inherit from the out_entry !?
# - class gtf

# - raise "soft" error if length of seq and qualSeq is different in fastq_rec 
# - add a log file to the parse methods to get information e.g about the 
#   number of read instances


import re

class align_entry():

    # set kimura to NA because the kimura distance is not stored in the
    # information line, it is stored after the alignment and is set
    # if the respective line is read
    def __init__(self, line):
        
        self.line = line.strip()
        self.line_split = self.line.split(" ")
        self.occurence = 1
        # Elements that are not on the complement strand do not have
        # the strand information in their information line
        if len(self.line_split) == 14:

            first_part = self.line_split[:8]
            second_part = ['+']
            third_part = self.line_split[8:]
            self.line_split = first_part + second_part + third_part

        else:

            self.line_split[8] = "-"

        self.swscore = self.line_split[0]
        self.div = self.line_split[1]
        self.del_bp = self.line_split[2]
        self.ins_bp = self.line_split[3]
        self.query_seq = self.line_split[4]
        self.chr_nr = self.query_seq.replace("chr", "")
        self.start = self.line_split[5]
        self.end = self.line_split[6]
        self.strand = self.line_split[8]
        self.repeat = self.line_split[9]
        self.entry_id = self.line_split[-1]
        self.kimura_distance = "NA"
               
        self.set_seq_id()
        
    # To get acess to each TE respectivly a unique id is necessary
    def set_seq_id(self):

        self.split_repeat_name()

        self.seq_id = ("{}|{}|{}|{}|{}|{}|{}|{}|{}".format(self.query_seq,
                                                self.strand,
                                                self.start,
                                                self.end,
                                                self.order,
                                                self.superorder,
                                                self.family,
                                                self.entry_id,
                                                self.kimura_distance))
  

    # parses the kimura distance and add them to the seq id
    def add_kimura(self, line):
        
        line = line.strip()
        kimura_distance = line.split("=")[1][1:]
        self.kimura_distance = kimura_distance
        self.set_seq_id()

    # the repeat name is usually combined of Family#order/superfamily
    # this method splits the repeat name to get access to the 
    # respective chategory
    def split_repeat_name(self):

        if re.search(r'\#', self.repeat):
            
            self.family = self.repeat.split('#')[0]
            
            if re.search(r'/', self.repeat):
            
                self.order = self.repeat.split('#')[1].split('/')[0]
                self.superorder = self.repeat.split('#')[1].split('/')[1]
            
            else:
            
                self.order = self.repeat.split('#')[1]
                self.superorder = "NA"
        else:

            self.family = "NA"
            self.order = "NA"
            self.superorder = "NA"

    # to create a typical line for a bed file    
    def translate_to_bed(self):

        separator = "\t"
            
        bed_line = separator.join((self.query_seq, 
                                   self.start,
                                   self.end,
                                   self.seq_id,
                                   self.swscore,
                                   self.strand))

        return bed_line


    def translate_to_gtf(self):

        separator = "\t"
        
        last_column = "\tgene_id \"{}\"; transcript_id \"{}\"; family_id \"{}\"; class_id \"{}\";".format(self.seq_id,
                        self.family,
                        self.order,
                        self.superorder)
   

        gtf_line = separator.join((self.chr_nr,
            "RepeatMasker",
             "exon",
             self.start,
             self.end,
             self.swscore,
             self.strand,
             last_column))

        return gtf_line


class bed_entry():

    def __init__(self, line):

        self.line = line.strip()
        self.line_split = self.line.split("\t")

        self.chrom = self.line_split[0]
        self.chromStart = int(self.line_split[1])
        self.chromEnd = int(self.line_split[2])
        self.length = self.chromEnd - self.chromStart

        # Iam not sure if the direction is considered 
        # so that I added this if call if the start 
        # coordinate is bigger than the end coordinate
        if self.length < 0:
            self.length = self.length * -1

        if len(self.line_split) > 3:
            self.name = self.line_split[3]
            self.score = self.line_split[4]
            self.strand = self.line_split[5]

    def __str__(self):
        
        separator = "\t"
        self.bed_line = separator.join((self.chrom, 
                                        str(self.chromStart),
                                        str(self.chromEnd),
                                        self.name,
                                        self.score,
                                        self.strand)) 
        
        return self.bed_line


    def write_to_bed_file(self, handle):

        separator = "\t"
        bed_line = separator.join((self.chrom, 
                                    str(self.chromStart),
                                    str(self.chromEnd),
                                    self.name,
                                    self.score,
                                    self.strand)) 

        handle.write(bed_line + "\n")


class seq_rec():

    def __init__(self, seqId, seq):
        self.seqId = seqId
        self.seq = seq


class fastq_rec(seq_rec):

    def __init__(self, seqId, seq, qualSeq):

        if len(seq) != len(qualSeq):
            print("sequence and quality string differ in length")
            seq_rec.__init__(self, "NA", "NA")
            self.qualSeq = "NA"
        else:
            seq_rec.__init__(self, seqId, seq)
            self.qualSeq = qualSeq

    def __str__(self):

        return self.seqId +"\n" \
                +self.seq + \
                "\n+\n" \
                +self.qualSeq


#
# Parse the information of each alignment and add the kimura distance
# 
def parse_align_file(file):

    align_entries_list = []

    with open(file, "r") as input_handle:

        for line in input_handle:
        
            if re.search(r'^[0-9]', line):

                te_instance = align_entry(line)

                align_entries_list.append(te_instance)
            
            if re.search(r'Kimura', line):
       
                te_instance.add_kimura(line)

                #yield te_instance
    
    return align_entries_list


def parse_bed_file(file):

    bed_entries_list = []

    with open(file, "r") as input_handle:

        for line in input_handle:

            annotation = bed_entry(line)

            bed_entries_list.append(annotation)

    return bed_entries_list


def parse_fasta_file(file):

    pass

# extracts all information from the align file and
# creats a bed line for each instance and stores them
# into a bed file
def align_to_bed(align_file, output_path):

    with open("%s/%s.bed" %(output_path, align_file.split(".")[0]), "w") as output_handle:
        
        for instance in parse_align_file(align_file):
           
            bed = bed_entry(instance.translate_to_bed()) 
           
            bed.write_to_bed_file(output_handle)


