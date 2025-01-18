import csv
import pandas as pd
from datetime import datetime
import os

def convert_to_csv(input_lines, filename):
    data = [line.strip().split('\t') for line in input_lines if line.strip()]

    headers = ['RT_len', 'RT_seq', 'RT_picked', 'PBS_len', 'PBS_seq', 'PBS_picked',
              '3_extension_seq', '3_extension_picked', 'extensF_oligo', 'extensR_oligo',
              'sgRNA_seq', 'sgRNA_rank', 'sgF_oligo', 'sgR_oligo', 'sg_Orientation',
              'sg_Seed/PAM_disrupt', 'sg_GC%', 'sg_OnTargetScore', 'Enzyme']

    with open(filename, 'w', newline='') as csvfile:
        writer = csv.writer(csvfile)
        writer.writerow(headers)
        writer.writerows(data)

def filter_pegrnas(input_csv, output_csv):
    df = pd.read_csv(input_csv)

    numeric_columns = ['PBS_len', 'RT_len', 'sgRNA_rank']
    for col in numeric_columns:
        df[col] = pd.to_numeric(df[col], errors='coerce')

    df = df.dropna(subset=numeric_columns)

    # Function to get filtered data for a specific rank
    def get_rank_filtered(df, rank):
        return df[
            (df['PBS_len'].between(12, 13)) &
            (df['RT_len'].between(13, 15)) &
            (~df['RT_seq'].str.startswith('C')) &
            (df['sgRNA_rank'] == rank) &
            (df['Enzyme'] == 'Cas9-NGG')
        ]

    # Process each rank
    results = {}
    for rank in [1, 2]:
        filtered_df = get_rank_filtered(df, rank)
        combinations = filtered_df.groupby(['PBS_len', 'RT_len']).size().reset_index(name='count')
        results[rank] = {
            'df': filtered_df,
            'combinations': combinations,
            'total': len(filtered_df)
        }

    # Print results with clear formatting
    print("\n" + "="*50)
    print("pegRNA Design Combinations Analysis")
    print("="*50)

    for rank in [1, 2]:
        print(f"\nRank {rank} sgRNAs:")
        print("-"*20)
        print(f"Total candidates: {results[rank]['total']}")
        print("\nPBS_len x RT_len combinations:")
        if len(results[rank]['combinations']) > 0:
            print(results[rank]['combinations'].to_string(index=False))
        else:
            print("No combinations found")
        print("-"*50)

    # Save combined filtered results
    combined_filtered = pd.concat([results[r]['df'] for r in [1, 2]])
    combined_filtered.to_csv(output_csv, index=False)
    # how to sort just a subsection of the rows? 
    # e.x. sort just the rows with sgRNA_rank = 0 by key1,
    # and sort rows with sgRNA_rank = 1 with key1.

    return results

import peglit
SCAFFOLD = 'GTTTTAGAGCTAGAAATAGCAAGTTAAAATAAGGCTAGTCCGTTATCAACTTGAAAAAGTGGCACCGAGTCGGTGC'

def get_linker(epegRNA):
    # spacer = sgRNA
    spacer, scaffold, template, pbs, motif = epegRNA
    linkers = peglit.pegLIT(
                seq_spacer=spacer,
                seq_scaffold=scaffold,
                seq_template=template,
                seq_pbs=pbs,
                # 3' motif
                seq_motif=motif
            )
    print(linkers)
    return linkers

def get_timestamped_filename(base_name):
    timestamp = datetime.utcnow().strftime('%Y%m%d_%H%M')
    return f"out/{timestamp}_{base_name}"

def ensure_out_directory():
    if not os.path.exists('out'):
        os.makedirs('out')

def main():
    ensure_out_directory()
    input_file = 'design.txt'
    with open(input_file, 'r') as file:
        text_content = file.read()
    
    lines = [line for line in text_content.split('\n') if line.strip()]
    
    initial_csv = get_timestamped_filename("pegRNAs.csv")
    convert_to_csv(lines, initial_csv)
    
    df = pd.read_csv(initial_csv)
    df = df.sort_values(by='sgRNA_rank')
    sorted_csv = get_timestamped_filename("sortedPeg.csv")
    df.to_csv(sorted_csv, index=False)
    
    filtered_csv = get_timestamped_filename("filtered_pegRNAs.csv")
    results = filter_pegrnas(sorted_csv, filtered_csv)
    
    # Read the filtered results
    filtered_df = pd.read_csv(filtered_csv)
    
    # Add new column for linkers
    filtered_df['linker'] = None
    
    # Process each row through get_linker
    for index, row in filtered_df.iterrows():
        epegRNA = (
            row['sgRNA_seq'],
            SCAFFOLD,
            row['RT_seq'],
            row['PBS_seq'],
            row['3_extension_seq']
        )
        linker = get_linker(epegRNA)[0]
        filtered_df.at[index, 'linker'] = linker
    
    # Save the results with timestamp
    final_output = get_timestamped_filename("pegRNAs_with_linkers.csv")
    filtered_df.to_csv(final_output, index=False)
    
    return results
if __name__ == "__main__":
    filtered_results = main()
