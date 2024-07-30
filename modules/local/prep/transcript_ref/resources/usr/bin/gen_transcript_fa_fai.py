import pysam
import pandas as pd
import subprocess

# Generate name dictionary for gencode to Verily chromosome/contig names from the Verily fa.fai file
def name_match(fai_file):
    names = set()
    with open(fai_file, 'r') as f:
        for line in f:
            name = line.split('\t')[0]
            names.add(name)
    return names

# Parse gft, 1. only keep 'exon' 2. generate a column with transcript_id 3. only keep entries where the chromosome/contig exits in Verily reference fasta
def parse_gtf(gtf_file, names):
    columns = [
        "seqname", "source", "feature", "start", "end", "score", "strand", "frame", "attribute"
    ]
    df = pd.read_csv(gtf_file, sep="\t", comment='#', header=None, names=columns)
    df = df[df['feature'] == 'exon']  # Filter for exon features only
    df['transcript_id'] = df['attribute'].str.extract('transcript_id "([^"]+)"')
    df = df[df['seqname'].isin(names)] # only keep the entries where the chromosome/contig exits in Verily
    return df

# Extract transcript sequences
def extract_transcript_sequences(fasta_file, gtf_df, names):
    fasta = pysam.FastaFile(fasta_file)
    transcripts = {}

    for _, row in gtf_df.iterrows():
        seqname = row['seqname']
        start = row['start'] - 1  # Convert to 0-based indexing
        end = row['end']
        transcript_id = row['transcript_id']

        if transcript_id not in transcripts:
            transcripts[transcript_id] = []

        exon_seq = fasta.fetch(seqname, start, end)
        transcripts[transcript_id].append(exon_seq)

    fasta.close()

    return transcripts

# Write transcripts in output file
def write_transcripts_to_fasta(transcripts, output_file):
    with open(output_file, 'w') as f:
        for transcript_id, exons in transcripts.items():
            full_sequence = ''.join(exons)
            f.write(f'>{transcript_id}\n')
            for i in range(0, len(full_sequence), 100):  # Format sequence to 60 characters per line
                f.write(full_sequence[i:i+100] + '\n')

def index_fasta(input_file):
    command = ["samtools", "faidx", input_file]
    subprocess.run(command)

def main():
    gtf_file = 'gencode.v40.chr_patch_hapl_scaff.annotation.gtf'
    fasta_file = 'GRCh38_Verily_v1.genome.fa'
    fai_file = 'GRCh38_Verily_v1.genome.fa.fai'
    output_file = 'GRCh38_Verily_v1.gencode_v40.transcriptome.exact00.fa'

    names = name_match(fai_file)
    gtf_df = parse_gtf(gtf_file, names)
    transcripts = extract_transcript_sequences(fasta_file, gtf_df, names)
    write_transcripts_to_fasta(transcripts, output_file)
    index_fasta(output_file)

if __name__ == "__main__":
    main()
