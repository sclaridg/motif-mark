#!/usr/bin/env python

# claridge_motif_mark.py

import re
import sys
import argparse
import cairo
import math
import random

parser = argparse.ArgumentParser(description="Creates a visualization of one to seven sequence motifs across a sequence or multiple sequences provided in UCSC format.\nRequires a FASTA file of the sequences to be visualized and a list of motifs to be visualized.\nSince a random color generator is used, if the output color scheme is unfavorable, run script again.")
parser.add_argument('-f','--fasta', help='absolute path to FASTA file of sequences to be searched (/path/to/<FASTA_FILE>).', required=True, type=str)
parser.add_argument('-m','--motifs', help='absolute path to the text file containing the list of motifs to be visualized, formatted as one motif per line (/path/to/<MOTIF_FILE>).', required=True, type=str)
args = parser.parse_args()


#####################################################################################
##### Define Higher Order Functions #################################################
#####################################################################################

##### Functions for file formatting and information

def fasta_twofer(filepath):
    '''Takes in a FASTA file and creates a new FASTA file called "sequences_twofer.fasta"
    that only has two lines for every FASTA entry.'''
    first_line = True
    with open(filepath,"r") as fasta, open("./sequences_twofer.fasta", "w+") as out:
        for line in fasta:
            if first_line and line.startswith(">"):     ### If the first line, print line without "\n" because print statement appends "\n"
                out.write(line.strip())
                out.write("\n")
                first_line = False     ### Never run through this part of the looping system again
            elif not first_line and line.startswith(">"):     ### Print "\n" preceeding all deflines
                out.write("\n")
                out.write(line.strip())
                out.write("\n")
            else:
                out.write(line.strip("\n"))     ### Strip "\n" from all sequence lines

def longest_sequence(filepath):
    '''Takes in a FASTA file and returns the length of the longest sequence in the file.'''
    with open(filepath,"r") as fasta:
        length = 0
        for line in fasta:
            line = line.strip()
            if line.startswith(">") == True:
                pass
            elif line.startswith(">") == False:
                new_length = len(line)
                if new_length > length:
                    length = new_length
                else:
                    pass
    return length
    
def count_lines(infile):
    '''Opens the input file and returns the number of lines in the file.'''
    with open(infile) as file:
        for i, line in enumerate(file):
            pass
    return i + 1

##### Functions for finding the motif positions

def motif2regex(filepath):
    '''Takes the list of motifs and returns a dictionary of regular expressions that will search for that motif in
    a given sequence, regardless of case. Keys are motifs; values are regex motifs.'''
    motif_dict = {}
    with open(filepath) as motif_list:
        for motif in motif_list:
            motif = motif.strip("\n")
            motif_upper = motif.upper()
            motif_components = list(motif_upper)
            for i in range(0, len(motif_components)):
                motif_components[i] = iupac.get(motif_components[i])
                regex_motif = "".join(motif_components)
                motif_dict[motif] = regex_motif
    return motif_dict

def lengths_intex(sequence):
    '''From a given UCSC sequence, i.e. lower case intron, upper case exon, lower case intron,
    this function returns a list of length four: length of the sequence, length of the first intron,
    length of the exon, and length of the second intron.'''
    sequence = sequence.strip()
    pattern = re.compile(r'([a-z]+)([A-Z]+)([a-z]+)')
    res = re.match(pattern, sequence)
    if len(res.groups()) != 3:
        raise ValueError("ERROR: Exiting program. Sequence is not in UCSC format.")
    lengths = [sequence, len(sequence), len(res.group(1)), len(res.group(2)), len(res.group(3))]
    return lengths

def sequence_parser(filepath):
    with open(filepath, "r+") as file:
        seq_details = {}
        header = ""
        sequence = ""
        first_line = True
        for line in file:
            if line.startswith(">") == True:
                line = line.strip()
                if first_line == False:
                    seq_details[header] = lengths_intex(sequence)
                header = line
                sequence = ""
            else:
                sequence = line
                first_line = False
        seq_details[header] = lengths_intex(sequence)
    return seq_details

def motif_position(motif_dict, header, sequence):
    '''Takes a motif in regular-expression format and searches for all instances of that motif
    in the provided sequence. Returns a list of 2-item lists containing start positions and lengths.
    Since Python idexes starting at 0, 1 is added to the start position to account for pixels when making the image.'''
    positions = {}
    for motif in motif_dict.keys():   
        for match in re.finditer(motif_dict[motif], sequence):
            s = match.start()
            e = match.end()
            if header + "_" + motif in positions:
                positions[header + "_" + motif].append([s + 1, e - s])
            elif header + "_" + motif not in positions:
                positions[header + "_" + motif] = [[s + 1, e - s]]
    return positions
        
def get_all_positions(fasta_f, motif_f):
    seq_details = sequence_parser(fasta_f)
    motif_dict = motif2regex(motif_f)
    all_positions = []
    for info in seq_details.items():
        # info[0] is the header
        # info[1] is the list of details
        # info[1][0] is the sequence
        # info[1][1] is the sequence length
        # info[1][2] is the length of intron 1
        # info[1][3] is the length of the exon
        # info[1][4] is the length of intron 2
        positions = motif_position(motif_dict, header = info[0], sequence = info[1][0])
        positions["intron1"] = info[1][2]
        positions["exon"] = info[1][3]
        positions["intron2"] = info[1][4]
        all_positions.append(positions)
    return all_positions

##### Functions for drawing the introns, exons, and motifs

def draw_intex(intron1, exon, intron2, start):
    '''Draws the introns and exon of a sequence with one pixel corresponding to one base.'''
    # draw whole length of intron1 + exon + intron2
    context.set_source_rgb(0, 0, 0)
    context.set_line_width(1)
    context.move_to(start[0], start[1])
    context.line_to(start[0] + intron1 + exon + intron2, start[1])
    context.stroke()
    
    # draw wider exon
    context.set_line_width(14)
    context.move_to(start[0] + intron1, start[1])
    context.line_to(start[0] + intron1 + exon, start[1])
    context.stroke()
    
def draw_title(start, gene, info):
    '''Draws in the gene name and sequence information.'''
    # draw gene name
    context.set_font_size(15)    
    context.set_source_rgb(0, 0, 0)
    context.move_to(start[0] - 325, start[1] - 5)
    context.show_text(gene)
    
    # draw sequence information
    context.set_font_size(10)
    context.move_to(start[0] - 325, start[1] + 5)
    context.show_text(info)

def draw_legend(start, motif, r, g, b):
    '''Draws in the colored motif legend.'''
    context.set_source_rgb(r, g, b)
    context.set_line_width(18)
    context.move_to(start[0] - 325, start[1] + spacing)
    context.show_text(motif)

def draw_motif(start, motif_position):
    '''Draws in a motif based on the start position and length of the motif.
    One pixel corresponds to one base.'''
    # draw motif
    context.move_to(start[0] + motif_position[0], start[1])
    context.line_to(start[0] + motif_position[0] + motif_position[1], start[1])
    context.stroke()
    
    # draw indicator tick mark
    context.move_to(start[0] + motif_position[0] + 0.25 * motif_position[1], start[1] - 15)
    context.show_text("|")
    
#####################################################################################
##### Check and Organize Files ######################################################
#####################################################################################

# set limit for number of motifs; too many motifs gets really messy
if count_lines(args.motifs) > 7:
    raise ValueError("ERROR: Exiting program. Too many motifs were supplied (maximum 7).")

# convert input fasta to have two lines per entry
fasta_twofer(args.fasta)   

#####################################################################################
##### Make Dictionaries #############################################################
#####################################################################################

### IUPAC ambiguity codes and associated regular expressions for ignoring case
iupac = {"A":"[Aa]", "T":"[TtUu]", "C":"[Cc]", "G":"[Gg]", "U":"[TtUu]",
         "M":"[AaCc]", "R":"[AaGg]", "W":"[AaTtUu]",
         "S":"[CcGg]", "Y":"[CcTtUu]", "K":"[GgTtUu]",
         "V":"[AaCcGg]", "H":"[AaCcTtUu]", "D":"[AaGgTtUu]", "B":"[CcGgTtUu]",
         "N":"[AaTtCcGgUu]"}

### motif color dictionary
color_dict = {}
with open(args.motifs) as motif_list:
    for motif in motif_list:
        motif = motif.strip("\n")
        r = random.random()
        g = random.random()
        b = random.random()
        color_dict[motif] = [r, g, b]

####################################################################################
##### Find Motifs in Sequences ######################################################
#####################################################################################

all_positions = get_all_positions(fasta_f = "./sequences_twofer.fasta", motif_f = args.motifs)

#####################################################################################
##### Generate Visualizations #######################################################
#####################################################################################

# set up Cairo surface
surface_width = longest_sequence("./sequences_twofer.fasta") + 350 + 25
surface_height = count_lines("./sequences_twofer.fasta") / 2 * 200

surface = cairo.SVGSurface("./motif_mark.svg", surface_width, surface_height)
context = cairo.Context(surface)
context.set_line_width(1)
start = [350,100] # [x, y]

for one_position in all_positions: # for each sequence...
    # draw the introns and exon
    draw_intex(one_position["intron1"], one_position["exon"], one_position["intron2"], start)
    # initialize spacing variable for drawing the legend
    spacing = 20
    
    for motif_group in one_position.items(): # for each motif found in this sequence...
        # motif_group[0] is the header line with the attached motif
        # motif_group[1] is the list of motif starting positions and lengths
        
        if motif_group[0].startswith(">"):
            motif = motif_group[0].split("_")[1]
            details = motif_group[0].split("_")[0][1:]
            gene = details.split(" ")[0]
            info = " ".join(details.split(" ")[1:])
            
            # draw in the gene name and legend
            draw_title(start, gene, info)
            draw_legend(start, motif, color_dict[motif][0], color_dict[motif][1], color_dict[motif][2]) 
            
            # add spacing between legend entries
            spacing += 15
            
            for m in motif_group[1]:
                draw_motif(start, m)
                
    # add spacing before next sequence            
    start = [start[0], start[1] + 150]

surface.finish()

print("Finished.")
sys.exit()
