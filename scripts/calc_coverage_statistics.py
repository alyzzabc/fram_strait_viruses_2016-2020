import pandas as pd
import argparse

parser = argparse.ArgumentParser()
parser.add_argument("-i","--in_bedgraph", type=str, required=True)
parser.add_argument("-o","--out_coverage_statistics", type=str, required=True)

args = parser.parse_args()

def calculate_coverage_statistics(df):
    # Group by sequence and calculate the proportion of bases covered
    sequence_coverage = df.groupby(df[0])[2].apply(lambda x: (x >= 1).sum() / len(x)).reset_index()
    sequence_coverage.columns = ['Sequence', 'Proportion']

    # Calculate mean coverage for each sequence
    mean_coverage = df.groupby(df[0])[2].mean().reset_index()
    mean_coverage.columns = ['Sequence', 'Mean_Coverage']
    
    #sequence_lengths = df.groupby(df[0]).apply(lambda x: pd.Series([x[0].iloc[0], max(x[1]) - min(x[1]) + 1], index=['Sequence', 'Sequence_Length'])).reset_index(drop=True)

    sequence_lengths = df.groupby(df[0])[1].apply(lambda x: (max(x) - min(x) + 1)).reset_index()
    sequence_lengths.columns = ['Sequence', 'Sequence_Length']

    # Merge the dataframes on 'Sequence' column
    result = pd.merge(sequence_lengths, sequence_coverage, on='Sequence')
    result = pd.merge(result, mean_coverage, on='Sequence', how='outer')

    return result

# Example usage
bedgraph_file = pd.read_csv(args.in_bedgraph, sep='\t', header=None)
coverage_statistics=calculate_coverage_statistics(bedgraph_file)
coverage_statistics.to_csv(args.out_coverage_statistics, sep='\t', index=False)
