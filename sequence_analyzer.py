import argparse
import csv
from pathlib import Path


DNA_BASES = {"A", "T", "C", "G"}

CODON_TABLE = {
    "TTT": "F", "TTC": "F", "TTA": "L", "TTG": "L",
    "CTT": "L", "CTC": "L", "CTA": "L", "CTG": "L",
    "ATT": "I", "ATC": "I", "ATA": "I", "ATG": "M",
    "GTT": "V", "GTC": "V", "GTA": "V", "GTG": "V",
    "TCT": "S", "TCC": "S", "TCA": "S", "TCG": "S",
    "CCT": "P", "CCC": "P", "CCA": "P", "CCG": "P",
    "ACT": "T", "ACC": "T", "ACA": "T", "ACG": "T",
    "GCT": "A", "GCC": "A", "GCA": "A", "GCG": "A",
    "TAT": "Y", "TAC": "Y", "TAA": "*", "TAG": "*",
    "CAT": "H", "CAC": "H", "CAA": "Q", "CAG": "Q",
    "AAT": "N", "AAC": "N", "AAA": "K", "AAG": "K",
    "GAT": "D", "GAC": "D", "GAA": "E", "GAG": "E",
    "TGT": "C", "TGC": "C", "TGA": "*", "TGG": "W",
    "CGT": "R", "CGC": "R", "CGA": "R", "CGG": "R",
    "AGT": "S", "AGC": "S", "AGA": "R", "AGG": "R",
    "GGT": "G", "GGC": "G", "GGA": "G", "GGG": "G",
}


def parse_fasta(file_path):
    sequences = []
    current_name = None
    current_sequence = []

    with open(file_path, "r", encoding="utf-8") as file:
        for line in file:
            line = line.strip()

            if not line:
                continue

            if line.startswith(">"):
                if current_name:
                    sequences.append({
                        "name": current_name,
                        "sequence": clean_sequence("".join(current_sequence))
                    })

                current_name = line[1:].strip()
                current_sequence = []
            else:
                current_sequence.append(line)

    if current_name:
        sequences.append({
            "name": current_name,
            "sequence": clean_sequence("".join(current_sequence))
        })

    return sequences


def clean_sequence(sequence):
    return sequence.upper().replace(" ", "").replace("\n", "").replace("\r", "")


def validate_dna(sequence):
    invalid_bases = sorted(set(sequence) - DNA_BASES)
    return invalid_bases


def nucleotide_counts(sequence):
    return {
        "A": sequence.count("A"),
        "T": sequence.count("T"),
        "C": sequence.count("C"),
        "G": sequence.count("G")
    }


def gc_content(sequence):
    if not sequence:
        return 0

    gc_count = sequence.count("G") + sequence.count("C")
    return (gc_count / len(sequence)) * 100


def transcribe_dna_to_rna(sequence):
    return sequence.replace("T", "U")


def reverse_complement(sequence):
    complement_map = str.maketrans("ATCG", "TAGC")
    return sequence.translate(complement_map)[::-1]


def translate_dna(sequence):
    protein = []

    for i in range(0, len(sequence) - 2, 3):
        codon = sequence[i:i + 3]
        amino_acid = CODON_TABLE.get(codon, "?")
        protein.append(amino_acid)

    return "".join(protein)


def find_open_reading_frames(sequence):
    open_reading_frames = []

    for frame in range(3):
        protein = translate_dna(sequence[frame:])
        start_index = None

        for index, amino_acid in enumerate(protein):
            if amino_acid == "M" and start_index is None:
                start_index = index

            if amino_acid == "*" and start_index is not None:
                open_reading_frames.append({
                    "frame": frame + 1,
                    "start_codon_position": frame + start_index * 3 + 1,
                    "stop_codon_position": frame + index * 3 + 1,
                    "protein_length": index - start_index + 1,
                    "protein_sequence": protein[start_index:index + 1]
                })
                start_index = None

    return open_reading_frames


def analyze_sequence(record):
    name = record["name"]
    sequence = record["sequence"]
    invalid_bases = validate_dna(sequence)

    if invalid_bases:
        return {
            "name": name,
            "valid": False,
            "invalid_bases": ",".join(invalid_bases),
            "length": len(sequence),
            "gc_content": None,
            "a_count": None,
            "t_count": None,
            "c_count": None,
            "g_count": None,
            "rna_sequence": None,
            "reverse_complement": None,
            "protein_sequence": None,
            "orf_count": None
        }

    counts = nucleotide_counts(sequence)
    orfs = find_open_reading_frames(sequence)

    return {
        "name": name,
        "valid": True,
        "invalid_bases": "",
        "length": len(sequence),
        "gc_content": round(gc_content(sequence), 2),
        "a_count": counts["A"],
        "t_count": counts["T"],
        "c_count": counts["C"],
        "g_count": counts["G"],
        "rna_sequence": transcribe_dna_to_rna(sequence),
        "reverse_complement": reverse_complement(sequence),
        "protein_sequence": translate_dna(sequence),
        "orf_count": len(orfs)
    }


def write_summary_csv(results, output_path):
    fieldnames = [
        "name",
        "valid",
        "invalid_bases",
        "length",
        "gc_content",
        "a_count",
        "t_count",
        "c_count",
        "g_count",
        "orf_count"
    ]

    with open(output_path, "w", newline="", encoding="utf-8") as file:
        writer = csv.DictWriter(file, fieldnames=fieldnames)
        writer.writeheader()

        for result in results:
            writer.writerow({key: result[key] for key in fieldnames})


def print_report(results):
    print("\nDNA Sequence Analysis Report")
    print("=" * 40)

    for result in results:
        print(f"\nSequence: {result['name']}")
        print("-" * 40)
        print(f"Valid DNA sequence: {result['valid']}")
        print(f"Length: {result['length']} bases")

        if not result["valid"]:
            print(f"Invalid bases: {result['invalid_bases']}")
            continue

        print(f"GC Content: {result['gc_content']}%")
        print(f"A: {result['a_count']}")
        print(f"T: {result['t_count']}")
        print(f"C: {result['c_count']}")
        print(f"G: {result['g_count']}")
        print(f"Open Reading Frames Found: {result['orf_count']}")
        print(f"Protein Translation: {result['protein_sequence'][:80]}")

        if len(result["protein_sequence"]) > 80:
            print("Protein Translation shown partially because sequence is long.")


def main():
    parser = argparse.ArgumentParser(
        description="Analyze DNA sequences from a FASTA file."
    )

    parser.add_argument(
        "input",
        help="Path to the input FASTA file."
    )

    parser.add_argument(
        "--output",
        default="sequence_summary.csv",
        help="Path for the output CSV summary file."
    )

    args = parser.parse_args()

    input_path = Path(args.input)

    if not input_path.exists():
        print(f"Error: Could not find input file: {input_path}")
        return

    records = parse_fasta(input_path)

    if not records:
        print("Error: No FASTA records found.")
        return

    results = [analyze_sequence(record) for record in records]

    print_report(results)
    write_summary_csv(results, args.output)

    print(f"\nSummary CSV saved to: {args.output}")


if __name__ == "__main__":
    main()
