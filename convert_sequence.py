# Define the mapping of old values to new values
conversion_map = {
    ">ref|NC_001133| [org=Saccharomyces cerevisiae] [strain=S288C] [moltype=genomic] [chromosome=I]": ">chr1",
    ">CHR2.19980913 Chromosome II Sequence": ">chr2",
    ">CHR3.19980521 Chromosome III Sequence": ">chr3",
    ">CHR4.19990210 Chromosome IV Sequence": ">chr4",
    ">CHR5.19970727 Chromosome V Sequence": ">chr5",
    ">CHR6.19970727 Chromosome VI Sequence": ">chr6",
    ">CHR7.19970703 Chromosome VII Sequence": ">chr7",
    ">CHR8.19970727 Chromosome VIII Sequence": ">chr8",
    ">ref|NC_001141| [org=Saccharomyces cerevisiae] [strain=S288C] [moltype=genomic] [chromosome=IX]": ">chr9",
    ">CHR10.19970727 Chromosome X Sequence": ">chr10",
    ">CHR11.19970727 Chromosome XI Sequence": ">chr11",
    ">CHR12.19970730 Chromosome XII Sequence": ">chr12",
    ">ref|NC_001145| [org=Saccharomyces cerevisiae] [strain=S288C] [moltype=genomic] [chromosome=XIII]": ">chr13",
    ">ref|NC_001146| [org=Saccharomyces cerevisiae] [strain=S288C] [moltype=genomic] [chromosome=XIV]": ">chr14",
    ">CHR15.19970811 Chromosome XV Sequence": ">chr15",
    ">ref|NC_001148| [org=Saccharomyces cerevisiae] [strain=S288C] [moltype=genomic] [chromosome=XVI]": ">chr16",
    ">ref|NC_001224| [org=Saccharomyces cerevisiae] [strain=S288C] [moltype=genomic] [location=mitochondrion] [top=circular]": ">chrM"
}

# Input and output file names
input_file = 'S288C_reference_sequence_R9-1-1_19990210.fsa'
output_file = 'S288C_sequence_iRibo.fsa'

# Read the input file, process lines, and write to the output file
with open(input_file, 'r') as infile, open(output_file, 'w') as outfile:
    for line in infile:
        # Check if the whole line matches any key in the conversion map
        if line.strip() in conversion_map:
            # Write the mapped value to the output file
            outfile.write(conversion_map[line.strip()] + '\n')
        else:
            # Write the original line to the output file if it doesn't match any key
            outfile.write(line)

print("Conversion and output completed successfully.")
