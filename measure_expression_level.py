import pandas as pd
import argparse

# TODO: Make this flexible by adding options for other inputs that are not featureCounts such as from stringtie
parser = argparse.ArgumentParser(description='Compute TPM (transcripts per million) from featureCounts output')
parser.add_argument('--featureCounts_file',required=True,help='Input the path to the featureCounts output')
parser.add_argument('--outfile',required=True,help='Input the path to the output file after computing TPM')
parser.add_argument('--normalization',required=True,help='Input the types of normalization. Currently only support either tpm or fpkm')

args = parser.parse_args()

# Read in the featureCounts file and skip the first row
data = pd.read_csv(args.featureCounts_file, skiprows=[0], sep='\t')
# Rename the column so that the last column is called Count
data.columns = ['Geneid', 'Chr', 'Start', 'End', 'Strand', 'Length', 'Count']

if args.normalization == 'tpm':
    # Divide the read counts by the length of each gene in kilobases. This gives you reads per kilobase (RPK).
    data['RPK'] = data['Count']/(data['Length']/1000)
    # print (data.head())

    # Count up all the RPK values in a sample and divide this number by 1,000,000. This is your “per million” scaling factor.
    # print (data['RPK'].sum())
    per_million_scaling_factor = data['RPK'].sum()/1000000

    # Divide the RPK values by the “per million” scaling factor. This gives you TPM.
    data['TPM'] = data['RPK']/per_million_scaling_factor

    # print (data.head())

    # Save to a file
    data.to_csv(args.outfile, sep='\t', index=False)

elif args.normalization == 'fpkm':
    # Count up the total reads in a sample and divide that number by 1,000,000 – this is our “per million” scaling factor.
    per_million_scaling_factor = data['Count'].sum()/1000000
    # Divide the read counts by the “per million” scaling factor. This normalizes for sequencing depth, giving you reads per million (RPM)
    data['RPM'] = data['Count']/per_million_scaling_factor
    # Divide the RPM values by the length of the gene, in kilobases. This gives you RPKM.
    data['FPKM'] = data['RPM']/(data['Length']/1000)
    # Save to a file
    data.to_csv(args.outfile, sep='\t', index=False)
