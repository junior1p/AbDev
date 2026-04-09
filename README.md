# AbDev — Antibody Developability Assessment

> In-silico developability profiling for therapeutic antibodies and nanobodies.
> No GPU. No API keys. Python 3.9+. Results in under 60 seconds.

[![PyPI version](https://img.shields.io/pypi/v/abdev?color=blue)](https://pypi.org/project/abdev/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

---

## What is AbDev?

AbDev is a fully automated antibody developability assessment pipeline that produces an industry-standard scorecard from a single amino acid sequence. It covers three assessment layers:

| Layer | What it checks | Method |
|-------|---------------|--------|
| **1 — Liability Scanning** | Deamidation, isomerization, oxidation, glycosylation, unpaired cysteines, RGD motifs | Regex motif scanning, stratified by CDR vs. framework |
| **2 — TAP Profiling** | CDR length, hydrophobicity, charge asymmetry vs. 242 clinical therapeutics | Five TAP metrics (Raybould et al. *PNAS* 2019) |
| **3 — Thera-SAbDab Benchmark** | Nearest approved/clinical therapeutic antibody | k-mer sequence identity against Oxford OPIG database |

**Output:** A composite 0–100 developability score with traffic-light classification (🟢 GREEN / 🟡 AMBER / 🔴 RED) and actionable engineering recommendations.

---

## Quick Start

### Python API

```python
from abdev_pipeline import run_abdev

# Single nanobody (VHH)
result = run_abdev(
    sequence_input="EVQLVESGGGLVQPGGSLRLSCAASGFNIKDTYIHWVRQAPGKGLEWVARIYPTNGYTRYADSVKGRFTISADTSKNTAYLQMNSLRAEDTAVYYCSRWGGDGFYAMDYWGQGTLVTVSS",
    name="my_nanobody",
    out_dir="abdev_results",
)

# Full IgG (VH:VL)
result = run_abdev(
    sequence_input="EVQLVESGGGLVQPGGSLRLSCAAS...:DIQMTQSPSSLSASVGDRVTIT...",
    name="my_antibody",
    out_dir="abdev_results",
)
```

### CLI

```bash
# Single sequence
python abdev_pipeline.py --seq "EVQLVESGGGLVQPGGS..." --name my_vhh

# VH:VL pair
python abdev_pipeline.py --vhvl "VH_sequence:VL_sequence" --name my_igG

# Batch from FASTA
python abdev_pipeline.py --fasta antibodies.fasta

# Skip Thera-SAbDab download (offline mode)
python abdev_pipeline.py --seq "..." --skip-benchmark
```

### Web Interface

```bash
python abdev_api.py --port 5001
# Open http://localhost:5001
```

### Docker

```bash
docker build -t abdev .
docker run -p 5001:5001 abdev
```

---

## Installation

```bash
# 1. Install HMMER (required by abnumber/ANARCI for antibody numbering)
sudo apt-get install hmmer    # Ubuntu/Debian
# conda install -c bioconda hmmer -y   # or via conda

# 2. Install AbDev
pip install abdev   # Coming soon to PyPI

# Or install from source:
git clone https://github.com/junior1p/AbDev.git
cd AbDev
pip install -e . --break-system-packages

# 3. Verify
python -c "from abdev_pipeline import run_abdev; print('AbDev ready')"
```

**Mainland China mirror:**
```bash
pip install abnumber pandas numpy matplotlib requests scipy biopython \
  -i https://pypi.tuna.tsinghua.edu.cn/simple --break-system-packages -q
```

---

## Output

```
abdev_results/my_antibody/
├── scorecard_Chain1.png      # Visual scorecard (gauge + TAP bars + liability chart)
├── scorecard_Chain2.png
└── my_antibody_results.json  # Full machine-readable assessment
```

### Scorecard

The scorecard has three panels:
1. **Composite score gauge** — semicircle with needle, color-coded GREEN/AMBER/RED
2. **TAP metrics bars** — five physicochemical metrics vs. clinical-stage therapeutic (CST) reference
3. **Liability chart** — stacked bar chart of chemical liabilities by region and severity

### Developability Score

| Score | Classification | Action |
|-------|--------------|--------|
| **75–100** | 🟢 GREEN | Good developability. Advance to biophysical characterization. |
| **50–74** | 🟡 AMBER | Moderate risk. Address flagged liabilities before lead selection. |
| **0–49** | 🔴 RED | High risk. Engineering campaign required before advancing. |

---

## Developability Metrics Reference

| Metric | HIGH Risk Threshold | Source |
|--------|--------------------|--------|
| CDR total length | > 55 aa | Raybould et al. 2019 |
| CDR hydrophobicity | > 2.5 (KD sum) | Raybould et al. 2019 |
| CDR positive charge | > 5.0 | Raybould et al. 2019 |
| CDR negative charge | < −4.0 | Raybould et al. 2019 |
| Charge asymmetry VH-VL | > ±2.0 | Raybould et al. 2019 |
| Asn-Gly (NG) in CDR | Any | Lu et al. mAbs 2019 |
| Asp-Gly (DG) in CDR | Any | Industry standard |
| N-X-S/T glycosylation | Any in variable domain | Industry standard |
| Unpaired cysteine | Any | Industry standard |

---

## Project Structure

```
AbDev/
├── abdev_pipeline.py       # Core pipeline (all 7 steps)
├── abdev_api.py           # Flask REST API
├── templates/
│   └── index.html         # Web interface
├── SKILL.md               # Claude Code skill (full natural-language instructions)
├── requirements.txt
├── README.md
└── LICENSE
```

---

## Scientific References

1. Raybould, M.I.J. et al. (2019). Five computational developability guidelines for therapeutic antibody profiling. *PNAS* 116 (24).
2. Raybould, M.I.J. et al. (2020). Thera-SAbDab: the Therapeutic Structural Antibody Database. *Nucleic Acids Research* 48 (D1).
3. Lu, X. et al. (2019). Deamidation and isomerization liability analysis of 131 clinical-stage antibodies. *mAbs* 11 (5).
4. Dunbar, J. & Deane, C.M. (2016). ANARCI: antigen receptor numbering and receptor classification. *Bioinformatics* 32 (4).
5. Sharma, V.K. et al. (2023). Blueprint for antibody biologics developability. *mAbs* 15 (1).

---

## License

MIT License — free for academic and commercial use.
