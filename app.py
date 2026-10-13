from __future__ import annotations

import csv
import io
from typing import Any

from flask import Flask, jsonify, render_template, request, send_file

app = Flask(__name__)

CODON_TABLE = {
    "UUU":"F","UUC":"F","UUA":"L","UUG":"L","UCU":"S","UCC":"S","UCA":"S","UCG":"S",
    "UAU":"Y","UAC":"Y","UAA":"*","UAG":"*","UGU":"C","UGC":"C","UGA":"*","UGG":"W",
    "CUU":"L","CUC":"L","CUA":"L","CUG":"L","CCU":"P","CCC":"P","CCA":"P","CCG":"P",
    "CAU":"H","CAC":"H","CAA":"Q","CAG":"Q","CGU":"R","CGC":"R","CGA":"R","CGG":"R",
    "AUU":"I","AUC":"I","AUA":"I","AUG":"M","ACU":"T","ACC":"T","ACA":"T","ACG":"T",
    "AAU":"N","AAC":"N","AAA":"K","AAG":"K","AGU":"S","AGC":"S","AGA":"R","AGG":"R",
    "GUU":"V","GUC":"V","GUA":"V","GUG":"V","GCU":"A","GCC":"A","GCA":"A","GCG":"A",
    "GAU":"D","GAC":"D","GAA":"E","GAG":"E","GGU":"G","GGC":"G","GGA":"G","GGG":"G"
}
STOP_CODONS = {"UAA", "UAG", "UGA"}


def clean_dna(value: Any) -> str:
    if not isinstance(value, str):
        raise ValueError("Enter a DNA sequence as text.")
    lines = value.splitlines()
    headers = [line for line in lines if line.lstrip().startswith(">")]
    if len(headers) > 1:
        raise ValueError("Use one FASTA record at a time.")
    sequence_text = "".join(line for line in lines if not line.lstrip().startswith(">"))
    sequence = "".join(sequence_text.upper().split())
    if not sequence:
        raise ValueError("Enter a DNA sequence.")
    invalid = sorted(set(sequence) - set("ATCG"))
    if invalid:
        raise ValueError("DNA sequence contains invalid character(s): " + ", ".join(invalid))
    if len(sequence) > 200_000:
        raise ValueError("For this demo, sequences are limited to 200,000 bases.")
    return sequence


def translate_dna(sequence: str) -> str:
    rna = sequence.replace("T", "U")
    amino_acids = []
    for i in range(0, len(rna) - 2, 3):
        aa = CODON_TABLE.get(rna[i:i+3])
        if aa is None or aa == "*":
            break
        amino_acids.append(aa)
    return "".join(amino_acids)


def analyze(sequence: str) -> dict:
    sequence = clean_dna(sequence)
    length = len(sequence)
    counts = {base: sequence.count(base) for base in "ATCG"}
    percentages = {base: round(counts[base] / length * 100, 4) for base in "ATCG"}
    complement = sequence.translate(str.maketrans("ATCG", "TAGC"))
    rna = sequence.replace("T", "U")
    return {
        "sequence": sequence,
        "length": length,
        "counts": counts,
        "percentages": percentages,
        "gc_content": round((counts["G"] + counts["C"]) / length * 100, 4),
        "complement": complement,
        "reverse_complement": complement[::-1],
        "rna": rna,
        "protein": translate_dna(sequence),
    }


def find_orfs(sequence: str) -> list[dict]:
    """Find the first in-frame stop after each AUG in all six frames."""
    sequence = clean_dna(sequence)
    results = []
    strands = [("+", sequence), ("-", sequence.translate(str.maketrans("ATCG", "TAGC"))[::-1])]
    for strand, dna in strands:
        rna = dna.replace("T", "U")
        for frame in range(3):
            open_starts = []
            for codon_start in range(frame, len(rna) - 2, 3):
                codon = rna[codon_start:codon_start + 3]
                if codon == "AUG":
                    open_starts.append(codon_start)
                elif codon in STOP_CODONS:
                    for start in open_starts:
                        orf_dna = dna[start:codon_start + 3]
                        results.append({
                            "strand": strand,
                            "frame": frame + 1,
                            "start": start + 1,
                            "end": codon_start + 3,
                            "length": codon_start + 3 - start,
                            "sequence": orf_dna,
                            "protein": translate_dna(orf_dna),
                        })
                    open_starts.clear()
    results.sort(key=lambda item: (item["strand"], item["start"], item["frame"]))
    return results


@app.get("/")
def home():
    return render_template("index.html")


@app.get("/api/health")
def health():
    return jsonify({"ok": True, "app": "GeneScope"})


@app.post("/api/analyze")
def api_analyze():
    try:
        data = request.get_json(silent=True) or {}
        return jsonify(analyze(data.get("sequence", "")))
    except ValueError as exc:
        return jsonify({"error": str(exc)}), 400


@app.post("/api/motif")
def api_motif():
    try:
        data = request.get_json(silent=True) or {}
        sequence = clean_dna(data.get("sequence", ""))
        motif = clean_dna(data.get("motif", ""))
        positions, start = [], 0
        while True:
            index = sequence.find(motif, start)
            if index < 0:
                break
            positions.append(index + 1)  # Human-friendly 1-based position
            start = index + 1  # Preserve overlapping hits
        return jsonify({"motif": motif, "positions": positions, "count": len(positions)})
    except ValueError as exc:
        return jsonify({"error": str(exc)}), 400


@app.post("/api/gc-window")
def api_gc_window():
    try:
        data = request.get_json(silent=True) or {}
        sequence = clean_dna(data.get("sequence", ""))
        window = int(data.get("window", 20))
        step = int(data.get("step", 1))
        if window < 1 or window > len(sequence):
            raise ValueError("Window size must be between 1 and the sequence length.")
        if step < 1:
            raise ValueError("Step size must be at least 1.")
        # Prefix counts make each window lookup constant time, even for a
        # large window and small step. This avoids recounting overlapping bases.
        gc_prefix = [0]
        for base in sequence:
            gc_prefix.append(gc_prefix[-1] + (base in "GC"))
        points = []
        for start in range(0, len(sequence) - window + 1, step):
            gc = (gc_prefix[start + window] - gc_prefix[start]) / window * 100
            points.append({"position": start + 1, "gc": round(gc, 4)})
        return jsonify({"window": window, "step": step, "points": points})
    except (ValueError, TypeError) as exc:
        return jsonify({"error": str(exc) or "Window and step must be whole numbers."}), 400


@app.post("/api/compare")
def api_compare():
    try:
        data = request.get_json(silent=True) or {}
        first = clean_dna(data.get("sequence1", ""))
        second = clean_dna(data.get("sequence2", ""))
        if len(first) != len(second):
            raise ValueError("For this simple comparison, both sequences must have the same length.")
        matches = [i for i, (a, b) in enumerate(zip(first, second), start=1) if a == b]
        mismatches = [i for i, (a, b) in enumerate(zip(first, second), start=1) if a != b]
        return jsonify({
            "length": len(first), "matches": len(matches), "mismatches": len(mismatches),
            "identity": round(len(matches) / len(first) * 100, 4),
            "match_positions": matches, "mismatch_positions": mismatches,
            "sequence1": first, "sequence2": second
        })
    except ValueError as exc:
        return jsonify({"error": str(exc)}), 400


@app.post("/api/orfs")
def api_orfs():
    try:
        data = request.get_json(silent=True) or {}
        return jsonify({"orfs": find_orfs(data.get("sequence", ""))})
    except ValueError as exc:
        return jsonify({"error": str(exc)}), 400


@app.post("/api/report")
def api_report():
    try:
        data = request.get_json(silent=True) or {}
        result = analyze(data.get("sequence", ""))
        output = io.StringIO()
        writer = csv.writer(output)
        writer.writerow(["Metric", "Value"])
        writer.writerow(["Sequence", result["sequence"]])
        writer.writerow(["Length (bases)", result["length"]])
        writer.writerow(["GC content (%)", result["gc_content"]])
        for base in "ATCG":
            writer.writerow([f"{base} count", result["counts"][base]])
            writer.writerow([f"{base} percentage", result["percentages"][base]])
        writer.writerow(["Complement", result["complement"]])
        writer.writerow(["Reverse complement", result["reverse_complement"]])
        writer.writerow(["RNA transcript", result["rna"]])
        writer.writerow(["Protein (translation from first base)", result["protein"]])
        memory = io.BytesIO(output.getvalue().encode("utf-8"))
        memory.seek(0)
        return send_file(memory, mimetype="text/csv", as_attachment=True, download_name="genescope_report.csv")
    except ValueError as exc:
        return jsonify({"error": str(exc)}), 400


if __name__ == "__main__":
    app.run(debug=True)
