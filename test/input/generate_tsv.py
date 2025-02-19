import sys
import csv
import os

def bedReadToIntronChain(line):
    """Extract intron chain from BED12 fields."""
    dir, start = line[5], int(line[1])
    esizes = [int(x) for x in line[10].split(',')[:-1]]
    estarts = [int(x) for x in line[11].split(',')[:-1]]

    introns = [
        (start + estarts[i] + esizes[i], start + estarts[i + 1])
        for i in range(len(esizes) - 1)
    ]
    return [x[::-1] for x in introns[::-1]] if dir == '-' else introns

def parse_bed_file(bed_path):
    """Parse BED file and return list of (intron_chain, start, end, full_line)."""
    transcript_data = []
    with open(bed_path, 'r') as bed:
        for line in bed:
            fields = line.strip().split('\t')
            intron_chain = bedReadToIntronChain(fields)
            start, end = int(fields[1]), int(fields[2])
            transcript_data.append((intron_chain, start, end, line.strip()))
    return transcript_data

def is_within_100bp(ref_start, ref_end, query_start, query_end):
    """Check if transcript ends are within 100 bp."""
    return abs(ref_start - query_start) <= 100 and abs(ref_end - query_end) <= 100

def transcript_matches(reference_transcripts, query_transcript):
    """Check if a query transcript matches reference transcripts in junctions and ends."""
    query_chain, query_start, query_end, _ = query_transcript
    return any(
        ref_chain == query_chain and is_within_100bp(ref_start, ref_end, query_start, query_end)
        for ref_chain, ref_start, ref_end, _ in reference_transcripts
    )

def generate_comparison_tsv(sample_folder, gt_bed, output_tsv):
    """Generate precision/recall metrics and write to a TSV file."""
    sample_name = os.path.basename(sample_folder.rstrip('/'))
    flair_bed = os.path.join(sample_folder, "collapsed.isoforms.bed")

    gt_transcripts = parse_bed_file(gt_bed)
    flair_transcripts = parse_bed_file(flair_bed)

    gt_introns = [chain for chain, _, _, _ in gt_transcripts]
    flair_introns = [chain for chain, _, _, _ in flair_transcripts]

    def compare_junction_chains(reference_chains, query_chains):
        matches = sum(1 for ref_chain in reference_chains if ref_chain in query_chains)
        return (matches / len(reference_chains)) * 100 if reference_chains else 0

    def compare_transcript_ends(reference_ends, query_ends):
        matches = sum(
            1 for ref_end in reference_ends
            if any(is_within_100bp(ref_end[0], ref_end[1], qe[0], qe[1]) for qe in query_ends)
        )
        return (matches / len(reference_ends)) * 100 if reference_ends else 0

    def compare_junction_chain_and_ends(reference_transcripts, query_transcripts):
        matches = sum(
            1 for ref_chain, ref_start, ref_end, _ in reference_transcripts
            if any(
                query_chain == ref_chain and is_within_100bp(ref_start, ref_end, query_start, query_end)
                for query_chain, query_start, query_end, _ in query_transcripts
            )
        )
        return (matches / len(reference_transcripts)) * 100 if reference_transcripts else 0

    junction_chain_recall = compare_junction_chains(gt_introns, flair_introns)
    junction_chain_precision = compare_junction_chains(flair_introns, gt_introns)

    gt_ends = [(start, end) for _, start, end, _ in gt_transcripts]
    flair_ends = [(start, end) for _, start, end, _ in flair_transcripts]

    transcript_ends_recall = compare_transcript_ends(gt_ends, flair_ends)
    transcript_ends_precision = compare_transcript_ends(flair_ends, gt_ends)

    junction_chain_and_ends_recall = compare_junction_chain_and_ends(gt_transcripts, flair_transcripts)
    junction_chain_and_ends_precision = compare_junction_chain_and_ends(flair_transcripts, gt_transcripts)

    with open(output_tsv, 'a', newline='') as tsvfile:
        writer = csv.writer(tsvfile, delimiter='\t')
        writer.writerow([
            sample_name,
            f"{junction_chain_recall:.2f}%",
            f"{junction_chain_precision:.2f}%",
            f"{transcript_ends_recall:.2f}%",
            f"{transcript_ends_precision:.2f}%",
            f"{junction_chain_and_ends_recall:.2f}%",
            f"{junction_chain_and_ends_precision:.2f}%"
        ])

def filter_transcripts(sample_folder, gt_bed, output_filtered_bed):
    """Filter FLAIR transcripts that match reference in both junction chains + ends."""
    flair_bed = os.path.join(sample_folder, "collapsed.isoforms.bed")
    if not os.path.exists(flair_bed):
        print(f"Error: {flair_bed} not found.")
        return

    gt_transcripts = parse_bed_file(gt_bed)
    flair_transcripts = parse_bed_file(flair_bed)

    matching_transcripts = [
        full_line for transcript in flair_transcripts
        if transcript_matches(gt_transcripts, transcript)
        for _, _, _, full_line in [transcript]
    ]

    with open(output_filtered_bed, 'w') as out_bed:
        for line in matching_transcripts:
            out_bed.write(line + '\n')

    print(f"Filtered BED file saved to: {output_filtered_bed} ({len(matching_transcripts)} transcripts retained)")

if __name__ == "__main__":
    if len(sys.argv) != 5:
        print("Usage: python script.py <sample_folder> <gt_bed> <output_tsv> <output_filtered_bed>")
        sys.exit(1)

    sample_folder = sys.argv[1]
    gt_bed = sys.argv[2]
    output_tsv = sys.argv[3]
    output_filtered_bed = sys.argv[4]

    if not os.path.exists(output_tsv):
        with open(output_tsv, 'w', newline='') as tsvfile:
            writer = csv.writer(tsvfile, delimiter='\t')
            writer.writerow([
                "samplename",
                "junction-chain-recall",
                "junction-chain-precision",
                "transcript-ends-recall",
                "transcript-ends-precision",
                "junction-chain-and-ends-recall",
                "junction-chain-and-ends-precision"
            ])

    generate_comparison_tsv(sample_folder, gt_bed, output_tsv)
    filter_transcripts(sample_folder, gt_bed, output_filtered_bed)
    print(f"Metrics written to {output_tsv} and filtered BED written to {output_filtered_bed}")
