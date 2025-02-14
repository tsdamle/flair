import pysam
import csv
import vcf
import os

def bedReadToIntronChain(line):
    dir, start, esizes, estarts = line[5], int(line[1]), [int(x) for x in line[10].split(',')[:-1]], [int(x) for x in line[11].split(',')[:-1]]
    introns = set()
    for i in range(len(esizes) - 1):
        introns.add((start + estarts[i] + esizes[i], start + estarts[i + 1]))
    if dir == '-':
        introns = {x[::-1] for x in introns}
    return introns

def parse_bed_file(bed_path):
    intron_chains = set()
    transcript_ends = set()
    with open(bed_path, 'r') as bed:
        for line in bed:
            fields = line.strip().split('\t')
            intron_chains.update(bedReadToIntronChain(fields))
            transcript_ends.add((int(fields[1]), int(fields[2])))
    return intron_chains, transcript_ends

def is_within_100bp(ref_start, ref_end, query_start, query_end):
    return abs(ref_start - query_start) <= 100 and abs(ref_end - query_end) <= 100

def compare_junction_chains(gt_chains, flair_chains):
    matches = len(gt_chains.intersection(flair_chains))
    return (matches / len(gt_chains)) * 100 if gt_chains else 0

def compare_transcript_ends(gt_ends, flair_ends):
    matches = sum(1 for fe in flair_ends if any(is_within_100bp(ge[0], ge[1], fe[0], fe[1]) for ge in gt_ends))
    return (matches / len(gt_ends)) * 100 if gt_ends else 0

def compare_junction_chain_and_ends(gt_transcripts, flair_transcripts):
    matches = 0
    for flair_chain, flair_start, flair_end in flair_transcripts:
        for gt_chain, gt_start, gt_end in gt_transcripts:
            if flair_chain == gt_chain and is_within_100bp(gt_start, gt_end, flair_start, flair_end):
                matches += 1
                break  # Ensure each flair transcript is counted only once
    return (matches / len(gt_transcripts)) * 100 if gt_transcripts else 0

def generate_comparison_tsv(sample_folder, gt_bed, output_tsv):
    sample_name = os.path.basename(sample_folder.rstrip('/'))
    flair_bed = os.path.join(sample_folder, "flair_collapse.bed")
    
    gt_introns, gt_ends = parse_bed_file(gt_bed)
    flair_introns, flair_ends = parse_bed_file(flair_bed)
    
    junction_chain_recall = compare_junction_chains(gt_introns, flair_introns)
    junction_chain_precision = compare_junction_chains(flair_introns, gt_introns)
    transcript_ends_recall = compare_transcript_ends(gt_ends, flair_ends)
    transcript_ends_precision = compare_transcript_ends(flair_ends, gt_ends)
    junction_chain_and_ends_recall = compare_junction_chain_and_ends(set(gt_introns), set(flair_introns))
    junction_chain_and_ends_precision = compare_junction_chain_and_ends(set(flair_introns), set(gt_introns))
    
    with open(output_tsv, 'a', newline='') as tsvfile:
        writer = csv.writer(tsvfile, delimiter='\t')
        writer.writerow([sample_name, f"{junction_chain_recall:.2f}%", f"{junction_chain_precision:.2f}%", f"{transcript_ends_recall:.2f}%", f"{transcript_ends_precision:.2f}%", f"{junction_chain_and_ends_recall:.2f}%", f"{junction_chain_and_ends_precision:.2f}%"])

# Write TSV header
with open("comparison.tsv", 'w', newline='') as tsvfile:
    writer = csv.writer(tsvfile, delimiter='\t')
    writer.writerow(["samplename", "junction-chain-recall", "junction-chain-precision", "transcript-ends-recall", "transcript-ends-precision", "junction-chain-and-ends-recall", "junction-chain-and-ends-precision"])
