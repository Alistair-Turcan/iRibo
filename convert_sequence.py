# Define the mapping of old values to new values
conversion_map = {
    ">ref|NC_001133| [org=Saccharomyces cerevisiae] [strain=S288C] [moltype=genomic] [chromosome=I]": ">chr1",
    ">ref|NC_001134| [org=Saccharomyces cerevisiae] [strain=S288C] [moltype=genomic] [chromosome=II]": ">chr2",
    ">ref|NC_001135| [org=Saccharomyces cerevisiae] [strain=S288C] [moltype=genomic] [chromosome=III]": ">chr3",
    ">ref|NC_001136| [org=Saccharomyces cerevisiae] [strain=S288C] [moltype=genomic] [chromosome=IV]": ">chr4",
    ">ref|NC_001137| [org=Saccharomyces cerevisiae] [strain=S288C] [moltype=genomic] [chromosome=V]": ">chr5",
    ">ref|NC_001138| [org=Saccharomyces cerevisiae] [strain=S288C] [moltype=genomic] [chromosome=VI]": ">chr6",
    ">ref|NC_001139| [org=Saccharomyces cerevisiae] [strain=S288C] [moltype=genomic] [chromosome=VII]": ">chr7",
    ">ref|NC_001140| [org=Saccharomyces cerevisiae] [strain=S288C] [moltype=genomic] [chromosome=VIII]": ">chr8",
    ">ref|NC_001141| [org=Saccharomyces cerevisiae] [strain=S288C] [moltype=genomic] [chromosome=IX]": ">chr9",
    ">ref|NC_001142| [org=Saccharomyces cerevisiae] [strain=S288C] [moltype=genomic] [chromosome=X]": ">chr10",
    ">ref|NC_001143| [org=Saccharomyces cerevisiae] [strain=S288C] [moltype=genomic] [chromosome=XI]": ">chr11",
    ">ref|NC_001144| [org=Saccharomyces cerevisiae] [strain=S288C] [moltype=genomic] [chromosome=XII]": ">chr12",
    ">ref|NC_001145| [org=Saccharomyces cerevisiae] [strain=S288C] [moltype=genomic] [chromosome=XIII]": ">chr13",
    ">ref|NC_001146| [org=Saccharomyces cerevisiae] [strain=S288C] [moltype=genomic] [chromosome=XIV]": ">chr14",
    ">ref|NC_001147| [org=Saccharomyces cerevisiae] [strain=S288C] [moltype=genomic] [chromosome=XV]": ">chr15",
    ">ref|NC_001148| [org=Saccharomyces cerevisiae] [strain=S288C] [moltype=genomic] [chromosome=XVI]": ">chr16",
    ">ref|NC_001224| [org=Saccharomyces cerevisiae] [strain=S288C] [moltype=genomic] [location=mitochondrion] [top=circular]": ">chrM"
}

# Input and output file names
input_file = 'S288C_reference_sequence_R64-3-1_20210421.fsa'
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
