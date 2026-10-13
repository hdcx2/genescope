"""Standalone DNA analyzer for the Project 1 command-line/Colab workflow.

Run with ``python dna_sequence_analyzer.py`` and enter a sequence, or use
``python dna_sequence_analyzer.py --self-test`` to check the sample cases.
The full interactive web app and its charts/tools are available through app.py.
"""

from __future__ import annotations

import sys

CODON_TABLE = {
    "UUU":"F", "UUC":"F", "UUA":"L", "UUG":"L", "UCU":"S", "UCC":"S", "UCA":"S", "UCG":"S",
    "UAU":"Y", "UAC":"Y", "UAA":"*", "UAG":"*", "UGU":"C", "UGC":"C", "UGA":"*", "UGG":"W",
    "CUU":"L", "CUC":"L", "CUA":"L", "CUG":"L", "CCU":"P", "CCC":"P", "CCA":"P", "CCG":"P",
    "CAU":"H", "CAC":"H", "CAA":"Q", "CAG":"Q", "CGU":"R", "CGC":"R", "CGA":"R", "CGG":"R",
    "AUU":"I", "AUC":"I", "AUA":"I", "AUG":"M", "ACU":"T", "ACC":"T", "ACA":"T", "ACG":"T",
    "AAU":"N", "AAC":"N", "AAA":"K", "AAG":"K", "AGU":"S", "AGC":"S", "AGA":"R", "AGG":"R",
    "GUU":"V", "GUC":"V", "GUA":"V", "GUG":"V", "GCU":"A", "GCC":"A", "GCA":"A", "GCG":"A",
    "GAU":"D", "GAC":"D", "GAA":"E", "GAG":"E", "GGU":"G", "GGC":"G", "GGA":"G", "GGG":"G",
}


def clean_dna(sequence: str) -> str:
    """Normalize case and whitespace, then reject empty or non-ATCG input."""
    lines = sequence.splitlines()
    headers = [line for line in lines if line.lstrip().startswith(">")]
    if len(headers) > 1:
        raise ValueError("Use one FASTA record at a time.")
    cleaned = "".join("".join(line for line in lines if not line.lstrip().startswith(">" )).upper().split())
    if not cleaned:
        raise ValueError("Please enter a DNA sequence.")
    invalid = sorted(set(cleaned) - set("ATCG"))
    if invalid:
        raise ValueError("DNA sequences may contain only A, T, C, and G. "
                         f"Invalid character(s): {', '.join(invalid)}")
    return cleaned


def analyze_dna(sequence: str) -> dict:
    """Return basic composition, strand conversions, and first-frame translation."""
    dna = clean_dna(sequence)
    length = len(dna)
    counts = {base: dna.count(base) for base in "ATCG"}
    percentages = {base: counts[base] / length * 100 for base in "ATCG"}
    complement = dna.translate(str.maketrans("ATCG", "TAGC"))
    rna = dna.replace("T", "U")
    protein = []
    for offset in range(0, len(rna) - 2, 3):
        amino_acid = CODON_TABLE[rna[offset:offset + 3]]
        if amino_acid == "*":
            break
        protein.append(amino_acid)
    return {
        "sequence": dna,
        "length": length,
        "counts": counts,
        "percentages": percentages,
        "gc_content": (counts["G"] + counts["C"]) / length * 100,
        "complement": complement,
        "reverse_complement": complement[::-1],
        "rna": rna,
        "protein": "".join(protein),
    }


def print_report(result: dict) -> None:
    """Display the analyzer results in a readable text report."""
    print("\nDNA SEQUENCE REPORT")
    print("=" * 40)
    print(f"Sequence:          {result['sequence']}")
    print(f"Length:            {result['length']} nucleotides")
    print(f"GC content:        {result['gc_content']:.2f}%")
    print("Nucleotide counts and percentages:")
    for base in "ATCG":
        print(f"  {base}: {result['counts'][base]:>6}  ({result['percentages'][base]:6.2f}%)")
    print(f"Complement:        {result['complement']}")
    print(f"Reverse complement:{result['reverse_complement']}")
    print(f"RNA transcript:    {result['rna']}")
    print(f"Protein:           {result['protein'] or '(no amino acids before a stop codon)'}")
    print("\nGC content is the percentage of bases that are guanine (G) or cytosine (C).")


def run_self_tests() -> None:
    """Check the four examples listed in the Project 1 instructions."""
    atgc = analyze_dna("ATGC")
    assert atgc["length"] == 4 and atgc["gc_content"] == 50
    assert atgc["counts"] == {"A": 1, "T": 1, "C": 1, "G": 1}
    assert analyze_dna("GGCC")["gc_content"] == 100
    try:
        analyze_dna("ATGX")
    except ValueError:
        pass
    else:
        raise AssertionError("ATGX should be rejected")
    stop_case = analyze_dna("ATGAAATAG")
    assert stop_case["protein"] == "MK"
    assert stop_case["counts"] == {"A": 5, "T": 2, "C": 0, "G": 2}
    assert round(stop_case["gc_content"], 2) == 22.22
    print("All four example checks passed.")


def main() -> None:
    if "--self-test" in sys.argv:
        run_self_tests()
        return
    try:
        print_report(analyze_dna(input("Enter a DNA sequence (A, T, C, G): ")))
    except ValueError as error:
        print(f"Error: {error}")
        raise SystemExit(1)


if __name__ == "__main__":
    main()
