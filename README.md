# GeneScope

GeneScope is an educational DNA sequence analyzer built with Python, Flask, HTML, CSS, and JavaScript. Open the page directly in a browser for client-side analysis, or run the Flask app locally.

## Features

- FASTA upload and pasted FASTA input (one record at a time)
- Base counts, percentages, sequence length, and GC content
- Complement, reverse complement, RNA transcript, and basic translation
- Nucleotide composition chart and sliding-window GC plot
- Motif search with overlapping matches
- Equal-length sequence comparison
- Candidate ORF search on both strands
- CSV report download
- Standalone Python command-line analyzer

## Open the website

No Python setup is needed for browser mode:

1. Open the `templates` folder.
2. Open `index.html` in your browser.
3. Paste DNA or upload one FASTA record, then click **Analyze sequence**.

The browser runs the analysis locally. The Flask backend is optional.

## Run with Flask

Install Python 3.10 or newer, then run these commands from the project folder.

### Windows

```powershell
py -m venv .venv
.venv\Scripts\activate
python -m pip install -r requirements.txt
python app.py
```

### macOS or Linux

```bash
python3 -m venv .venv
source .venv/bin/activate
python -m pip install -r requirements.txt
python app.py
```

Live Website: <a>https://genescope-9o6z.onrender.com</a>

## Command-line analyzer

The standalone script uses only the Python standard library:

```bash
python dna_sequence_analyzer.py
python dna_sequence_analyzer.py --self-test
```

## Example checks

- `ATGC` has length 4 and 50% GC content.
- `GGCC` has 100% GC content.
- `ATGX` is rejected as invalid input.
- `ATGAAATAG` translates to `MK` before the stop codon.
- Motif `ATA` occurs at positions 1 and 3 in `ATATAT` (overlapping matches count).

## Project files

```text
GeneScope/
├── app.py                    # Flask app and analysis endpoints
├── dna_sequence_analyzer.py  # Standalone command-line analyzer
├── requirements.txt          # Flask dependency
├── templates/
│   └── index.html            # Website markup
├── static/
│   ├── app.js                # Browser interactions and local analysis
│   └── style.css             # Responsive styles
├── .gitignore
└── README.md
```

## Notes and limits

This is an educational project, not a clinical tool. Translation starts at the first base and stops at the first stop codon; it does not determine the correct biological reading frame. ORFs are candidates, not confirmed genes. Sequence comparison is positional and requires equal-length inputs. FASTA files must contain a single record.

## Deploy as a live website

The repository includes `render.yaml` for deployment on Render. Push the project to GitHub, then in Render choose **New** → **Blueprint**, connect the repository, and let Render use the blueprint configuration. It installs the dependencies from `requirements.txt` and starts the Flask app with Gunicorn. Render will provide a public `onrender.com` URL when deployment finishes.

You can also create a Web Service manually with build command `pip install -r requirements.txt` and start command `gunicorn app:app`.
