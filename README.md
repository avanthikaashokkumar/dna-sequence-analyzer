# DNA Sequence Analyzer

A Python-based bioinformatics tool for analyzing DNA sequences from FASTA files.

## Project Overview

This project reads DNA sequences from a FASTA file and performs basic computational biology analysis. It validates DNA sequences, calculates nucleotide composition, estimates GC content, creates RNA transcripts, generates reverse complements, translates DNA into protein sequences, and identifies possible open reading frames.

This project was built to practice Python programming, bioinformatics logic, file parsing, and biological data analysis.

## Features

- Reads FASTA files
- Cleans and validates DNA sequences
- Detects invalid bases
- Calculates sequence length
- Counts A, T, C, and G nucleotides
- Calculates GC content
- Transcribes DNA to RNA
- Generates reverse complement sequences
- Translates DNA into protein sequences
- Detects possible open reading frames
- Exports summary results to CSV

## How to Run

Clone the repo and run:

```bash
python sequence_analyzer.py sample_sequences.fasta
```

To choose your own output file name:

```bash
python sequence_analyzer.py sample_sequences.fasta --output results.csv
```

## Example FASTA Format

```fasta
>sample_gene
ATGGCCATTGTAATGGGCCGCTGAAAGGGTGCCCGATAG
```

## Output

The program prints a report in the terminal and also creates a CSV file with summary results.

Example CSV columns:

```csv
name,valid,invalid_bases,length,gc_content,a_count,t_count,c_count,g_count,orf_count
```

## Biological Concepts Used

DNA is made of four nucleotide bases: adenine, thymine, cytosine, and guanine.

GC content measures the percentage of guanine and cytosine bases in a DNA sequence. GC content can affect DNA stability, gene structure, and sequencing behavior.

Transcription converts DNA into RNA by replacing thymine with uracil.

Translation converts DNA codons into amino acids using the genetic code.

Open reading frames are possible protein-coding regions that begin with a start codon and end with a stop codon.

## Why I Built This

I built this project to connect Python programming with bioinformatics fundamentals.

Many biological datasets are stored as sequence files, so learning how to parse, clean, analyze, and export sequence data is an important beginner skill for computational biology.

This project shows that I can use Python to process biological data, apply biological rules, and create reproducible outputs.

## How It Works

The program follows this workflow:

1. Reads DNA sequences from a FASTA file
2. Cleans each sequence by removing spaces and line breaks
3. Checks whether the sequence only contains valid DNA bases
4. Calculates nucleotide counts and GC content
5. Transcribes DNA into RNA
6. Generates the reverse complement sequence
7. Translates DNA codons into amino acids
8. Searches for possible open reading frames
9. Prints a terminal report
10. Saves a summary CSV file

## Files in This Repository

```text
dna-sequence-analyzer/
├── sequence_analyzer.py
├── sample_sequences.fasta
└── README.md
```

## Sample Dataset

The included `sample_sequences.fasta` file contains example DNA sequences that can be used to test the program.

One sequence intentionally contains an invalid base so the program can demonstrate error detection.

## Limitations

This is a beginner-friendly analysis tool and not a full professional genome annotation pipeline.

The translation function assumes a standard codon table and does not account for alternative genetic codes, introns, sequencing errors, or post-transcriptional modifications.

Open reading frame detection is simplified and should be interpreted as a learning feature rather than a validated gene prediction tool.

## Future Improvements

- Add FASTQ support
- Add sequence quality scoring
- Add visualization of nucleotide composition
- Add codon usage analysis
- Add ORF export to a separate file
- Add support for multiple genetic codes
- Add unit tests with pytest
- Add a simple Streamlit web app version

## Resume Bullet

Built a Python bioinformatics tool to parse FASTA files, validate DNA sequences, calculate GC content, translate DNA into protein sequences, detect possible open reading frames, and export results to CSV.

## GitHub Description

A Python bioinformatics tool for FASTA parsing, DNA validation, GC content analysis, transcription, translation, open reading frame detection, and CSV export.
