import os
from collections import defaultdict

# =========================================================
# FASTA UTILITIES
# =========================================================

def read_fasta(filepath):
    with open(filepath, 'r') as f:
        lines = f.readlines()

    sequence = ""
    for line in lines:
        if not line.startswith(">"):
            sequence += line.strip().upper()

    return sequence


def write_fasta(filepath, header, sequence):
    with open(filepath, 'w') as f:
        f.write(f">{header}\n")
        for i in range(0, len(sequence), 80):
            f.write(sequence[i:i+80] + "\n")


# =========================================================
# ORI DETECTION (Slide-Compliant)
# =========================================================

def skew_array(genome):
    skew = [0]
    for base in genome:
        if base == 'G':
            skew.append(skew[-1] + 1)
        elif base == 'C':
            skew.append(skew[-1] - 1)
        else:
            skew.append(skew[-1])
    return skew


def minimum_skew_positions(genome):
    skew = skew_array(genome)
    min_value = min(skew)
    return [i for i, v in enumerate(skew) if v == min_value]


def frequency_map(text, k):
    freq = defaultdict(int)
    for i in range(len(text) - k + 1):
        kmer = text[i:i+k]
        freq[kmer] += 1
    return freq


def reverse_complement(text):
    complement = {'A':'T','T':'A','C':'G','G':'C'}
    return ''.join(complement[b] for b in reversed(text))


def detect_ori(genome, window_size=500):
    print("\n=== ORI DETECTION ===")
    min_positions = minimum_skew_positions(genome)

    if not min_positions:
        raise Exception("Could not detect ORI using skew.")

    ori_center = min_positions[0]
    print(f"Predicted ORI position (minimum skew): {ori_center}")

    # Handle circular genome
    circular_genome = genome + genome
    ori_region = circular_genome[ori_center:ori_center + window_size]

    # Identify DnaA-like motifs (frequent 9-mers)
    freq = frequency_map(ori_region, 9)
    max_count = max(freq.values())
    frequent_kmers = [k for k, v in freq.items() if v == max_count]

    print("Most frequent 9-mers (possible DnaA boxes):")
    for kmer in frequent_kmers:
        print(f"{kmer} | RC: {reverse_complement(kmer)}")

    return ori_region


# =========================================================
# DESIGN + MARKER HANDLING
# =========================================================

def parse_design(filepath):
    components = []
    with open(filepath, 'r') as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("**"):
                continue
            key = line.split(',')[0].strip()
            components.append(key)
    return components


def parse_markers(filepath):
    markers = {}
    with open(filepath, 'r') as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            parts = line.split(',')
            if len(parts) >= 2:
                markers[parts[0].strip()] = parts[1].strip()
    return markers


# =========================================================
# MODE 1 — MODIFY EXISTING PLASMID
# =========================================================

def modify_existing_plasmid(sequence, design_components):
    print("\n=== MODIFYING EXISTING PLASMID ===")
    modified = sequence

    # EcoRI deletion logic
    if "EcoRI_site" not in design_components:
        modified = modified.replace("GAATTC", "")
        modified = modified.replace("CTTAAG", "")
        print("EcoRI site removed.")

    return modified


# =========================================================
# MODE 2 — BUILD NEW PLASMID FROM UNKNOWN GENOME
# =========================================================

def build_new_plasmid(genome, design_components, markers):
    print("\n=== BUILDING NEW PLASMID ===")

    ori_seq = detect_ori(genome)

    plasmid_seq = ori_seq

    for comp in design_components:
        if comp in markers:
            plasmid_seq += markers[comp]
            print(f"Added component: {comp}")
        else:
            print(f"Warning: {comp} not found in markers dictionary.")

    return plasmid_seq


# =========================================================
# MAIN
# =========================================================

def main():
    base_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    input_dir = os.path.join(base_dir, "Inputs")
    output_dir = os.path.join(base_dir, "Output")
    os.makedirs(output_dir, exist_ok=True)

    fasta_path = os.path.join(input_dir, "pUC19.fa")
    markers_path = os.path.join(input_dir, "markers.tab")
    design_path = os.path.join(input_dir, "Design_pUC19.txt")
    output_path = os.path.join(output_dir, "Output.Fa")

    sequence = read_fasta(fasta_path)
    design_components = parse_design(design_path)
    markers = parse_markers(markers_path)

    # Decide mode
    if "pUC19" in fasta_path:
        final_seq = modify_existing_plasmid(sequence, design_components)
    else:
        final_seq = build_new_plasmid(sequence, design_components, markers)

    write_fasta(output_path, "Modified_Plasmid", final_seq)

    print("\nFinal plasmid length:", len(final_seq))

    if "GAATTC" in final_seq or "CTTAAG" in final_seq:
        print("WARNING: EcoRI site still present.")
    else:
        print("EcoRI successfully removed.")

    print("Done.")


if __name__ == "__main__":
    main()
