# Define the mapping of old values to new values
conversion_map = {
    'chrI': 'chr1',
    'chrII': 'chr2',
    'chrIII': 'chr3',
    'chrIV': 'chr4',
    'chrIX': 'chr9',
    'chrmt': 'chrM',
    'chrV': 'chr5',
    'chrVI': 'chr6',
    'chrVII': 'chr7',
    'chrVIII': 'chr8',
    'chrX': 'chr10',
    'chrXI': 'chr11',
    'chrXII': 'chr12',
    'chrXIII': 'chr13',
    'chrXIV': 'chr14',
    'chrXV': 'chr15',
    'chrXVI': 'chr16'
}

# Input and output file names
input_file = 'saccharomyces_cerevisiae.gff'
output_file = 'saccharomyces_cerevisiae_iRibo.gff'

# Read the input file, process lines, and write to the output file
with open(input_file, 'r') as infile, open(output_file, 'w') as outfile:
    for line in infile:
        # Stop scanning the file when reaching the ##FASTA line
        if line.startswith('##FASTA'):
            break

        # Split the line by tab to access the first column
        columns = line.strip().split('\t')
        if len(columns) > 0:
            # Check if the first column value exists in the mapping, and replace it if it does
            if columns[0] in conversion_map:
                columns[0] = conversion_map[columns[0]]

        # Write the modified line to the output file
        outfile.write('\t'.join(columns) + '\n')

print("Conversion and output completed successfully.")
