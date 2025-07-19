import sys
from Bio import SeqIO
import pandas as pd

def calculate_content(sequence, window_size):
    results = []
    for i in range(0, len(sequence), window_size):
        window_seq = sequence[i:i + window_size]
        total_bases = len(window_seq)
        if total_bases == 0:
            continue

        lower_bases = sum(1 for base in window_seq if base.islower())
        gc_content = sum(1 for base in window_seq if base in 'GCgc')

        lower_ratio = lower_bases / total_bases
        gc_ratio = gc_content / total_bases

        results.append((i, lower_ratio, gc_ratio))

    return results

def process_genome(input_path, output_path, window_size):
    data = []

    for record in SeqIO.parse(input_path, "fasta"):
        chromosome = record.id
        sequence = str(record.seq)

        content = calculate_content(sequence, window_size)
        for start, lower_ratio, gc_ratio in content:
            data.append((chromosome, start, start + window_size, lower_ratio, gc_ratio))

    df = pd.DataFrame(data, columns=['Chromosome', 'Start', 'End', 'Lowercase Ratio', 'GC Ratio'])
    df.to_csv(output_path, index=False, sep='\t')
    print(f"Analysis complete. Results saved to {output_path}.")

if __name__ == "__main__":
    if len(sys.argv) != 4:
        print("Usage: python test.py input.fa output.txt window_size")
        sys.exit(1)

    input_file = sys.argv[1]
    output_file = sys.argv[2]
    window_size = int(sys.argv[3])

    process_genome(input_file, output_file, window_size)
