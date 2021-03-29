#!/gsc/biosw/bin/python3.7
import sys
import re

align_file = sys.argv[1]
out_file_entries = [] # the list is filled with strings that are stored in out

def delEntry(array, pos):
    
    return array[:pos] + array[pos+1:]
   
def writeOutFile(entries):

    with open("mm10.fa.out", "w") as output_handle:

        for entry in entries:
            output_handle.write(entry)

with open(align_file, "r") as input_handle:

    for line in input_handle:

        if re.search(r'^[0-9]', line):

           line = line.strip() 

           line = line.split(" ")

           line = delEntry(line, -2)
          
           if len(line) == 13:

               out_line = "{} {} {}".format(' '.join(line[:8]),
                       "+",
                       ' '.join(line[8:]))
           else:

               out_line = "{}".format(' '.join(line))

           out_line = ' '.join(out_line.split("#"))
           
          
           print(out_line)
            
