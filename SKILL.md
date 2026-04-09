# AbDev: Antibody Developability Assessment Skill

## Trigger
Use this skill when the user wants to:
- Screen an antibody or nanobody sequence for CMC (Chemistry, Manufacturing & Controls) risks
- Identify chemical liability hotspots: deamidation, isomerization, oxidation, glycosylation
- Compute physicochemical properties: CDR length, charge, hydrophobicity, pI, instability index
- Benchmark a candidate against FDA/EMA-approved therapeutic antibodies (Thera-SAbDab)
- Generate a developability scorecard for lead selection
- Compare multiple antibody variants side-by-side for developability triage

Example triggers:
- "Check this antibody sequence for developability issues"
- "Scan my nanobody VHH for deamidation and oxidation hotspots"
- "How does this antibody compare to approved therapeutics in terms of CDR length and charge?"
- "Give me a developability report for these 5 antibody variants"
- "Run AbDev on trastuzumab and flag any liabilities"

## Overview

**AbDev** is a fully automated antibody developability assessment pipeline that produces an industry-standard scorecard from a single amino acid sequence. It requires no GPU, no structural model, and no external API keys — only Python.

The pipeline covers **three assessment layers**:

**Layer 1 — Sequence Liability Scanning:**
Systematic motif search for chemical modification hotspots in the variable domain, stratified by CDR vs. framework region location (CDR liabilities are higher risk).

**Layer 2 — Physicochemical Profiling:**
Compute the five key metrics from the Therapeutic Antibody Profiler (TAP, Raybould et al. PNAS 2019): CDR length, surface hydrophobicity proxy, positive CDR charge, negative CDR charge, and charge asymmetry. Compare against the clinical-stage therapeutic (CST) reference distribution.

**Layer 3 — Thera-SAbDab Benchmarking:**
Query the Oxford OPIG Thera-SAbDab database (all WHO-recognised therapeutic antibodies) to find the most similar approved therapeutic, providing a "nearest approved neighbour" context for novelty and safety assessment.

**Output:** A one-page developability scorecard (PDF/text) with a 0–100 risk score, traffic-light classification (GREEN/AMBER/RED) per metric, and actionable engineering recommendations.

**Scientific rationale:**
Post-translational modifications such as deamidation, isomerization, and oxidation cause product heterogeneity and can significantly reduce antibody potency and stability. Early in-silico flagging of these liabilities is standard practice in industrial antibody discovery pipelines, enabling engineering fixes (e.g. Asn→Gln to prevent deamidation) before costly late-stage CMC failures.

---

## Step-by-Step Instructions for the Agent

### Step 0: Environment Setup

```bash
# Python 3.9+ required
python3 --version

# Install HMMER (system package — required by ANARCI/abnumber)
# Ubuntu/Debian:
sudo apt-get install -y hmmer

# Or via conda:
conda install -c bioconda hmmer=3.3.2 -y

# Install Python dependencies
pip install abnumber pandas numpy matplotlib requests scipy biopython \
  --break-system-packages -q

# Verify
python3 -c "from abnumber import Chain; print('AbNumber/ANARCI ready')"
python3 -c "import pandas, numpy, matplotlib, requests; print('Core libraries ready')"
```

Mainland China mirror:
```bash
pip install abnumber pandas numpy matplotlib requests scipy biopython \
  -i https://pypi.tuna.tsinghua.edu.cn/simple --break-system-packages -q
```

**Note on ANARCI/abnumber:** abnumber bundles ANARCI and requires HMMER to be installed at the system level. The `hmmer` package is available via `apt-get install hmmer` on Ubuntu/Debian, or via conda.

---

### Step 1: Input Parsing and Antibody Numbering

```python
from abnumber import Chain
import re

# ── IMGT CDR definitions (position ranges) ─────────────────────────────────
IMGT_CDR_POSITIONS = {
  "H": {"CDR1": range(27, 39), "CDR2": range(56, 66), "CDR3": range(105, 118)},
  "L": {"CDR1": range(27, 39), "CDR2": range(56, 66), "CDR3": range(105, 118)},
  "K": {"CDR1": range(27, 39), "CDR2": range(56, 66), "CDR3": range(105, 118)},
}

def parse_antibody_input(input_str):
  """
  Parse antibody sequence input.
  Accepts:
    - Single sequence string (VH or VHH/nanobody)
    - Two sequences separated by ':' (VH:VL)
    - FASTA format (single or two entries)
    - Named chains dict: {"VH": "EVQL...", "VL": "DIQM..."}
  Returns dict with chain annotations and IMGT numbering.
  """
  chains = {}
  if isinstance(input_str, dict):
    sequences = input_str
  elif ":" in input_str and not input_str.startswith(">"):
    parts = input_str.split(":", 1)
    sequences = {"Chain1": parts[0].strip(), "Chain2": parts[1].strip()}
  elif input_str.startswith(">"):
    sequences = {}
    current_name = None
    for line in input_str.strip().split("\n"):
      if line.startswith(">"):
        current_name = line[1:].strip().split()[0]
        sequences[current_name] = ""
      else:
        if current_name:
          sequences[current_name] += line.strip().upper()
  else:
    sequences = {"Chain1": input_str.strip().upper()}

  numbered_chains = {}
  for name, seq in sequences.items():
    seq = re.sub(r'[^ACDEFGHIKLMNPQRSTVWY]', '', seq.upper())
    if len(seq) < 50:
      print(f"Warning: {name} sequence too short ({len(seq)} aa)")
      continue
    try:
      chain = Chain(seq, scheme="imgt")
      chain_type = chain.chain_type  # "H", "L", "K" or "VHH"

      cdrs = {
        "CDR1": chain.cdr1_seq,
        "CDR2": chain.cdr2_seq,
        "CDR3": chain.cdr3_seq,
      }

      # Build framework sequence
      framework_parts = []
      for imgt_pos, aa in chain:
        in_cdr = any(
          imgt_pos.number in cdr_range
          for cdr_range in IMGT_CDR_POSITIONS.get(chain_type, {}).values()
        )
        if not in_cdr:
          framework_parts.append(aa)

      numbered_chains[name] = {
        "sequence": seq,
        "chain_type": chain_type,
        "chain_obj": chain,
        "cdrs": cdrs,
        "cdr_combined": "".join(cdrs.values()),
        "framework": "".join(framework_parts),
        "length": len(seq),
        "cdr_total_length": sum(len(v) for v in cdrs.values()),
      }

      ct_label = "VHH/Nanobody" if chain_type == "H" and len(seq) < 140 else f"V{chain_type}"
      print(f"{name}: {ct_label}, {len(seq)} aa, "
            f"CDR lengths: H1={len(cdrs.get('CDR1',''))}, "
            f"H2={len(cdrs.get('CDR2',''))}, "
            f"H3={len(cdrs.get('CDR3',''))}")

    except Exception as e:
      print(f"Warning: ANARCI numbering failed for {name}: {e}")
      print("Falling back to raw sequence analysis (no CDR annotation)")
      numbered_chains[name] = {
        "sequence": seq, "chain_type": "unknown",
        "chain_obj": None, "cdrs": {}, "cdr_combined": seq,
        "framework": seq, "length": len(seq), "cdr_total_length": 0,
      }

  return numbered_chains
```

---

### Step 2: Chemical Liability Scanning (Layer 1)

```python
import re
from dataclasses import dataclass

@dataclass
class Liability:
  name: str
  motif: str
  position: int       # 0-indexed in full sequence
  region: str         # "CDR1", "CDR2", "CDR3", "FR1-4", or "unknown"
  severity: str       # "HIGH" (CDR), "MEDIUM" (FR), "LOW"
  description: str
  recommendation: str

# ── Liability motif definitions ─────────────────────────────────────────────
# Source: industry practice + Lu et al. mAbs 2019, Raybould et al. PNAS 2019
LIABILITY_PATTERNS = [
  # Asn deamidation
  {"name": "Asn_Deamidation_NG",    "pattern": r"NG",
   "severity_cdr": "HIGH", "severity_fr": "MEDIUM",
   "description": "Asn-Gly: fastest deamidation motif (t½ ~1 day at pH 7.4)",
   "recommendation": "Mutate N→Q or N→A if in CDR paratope; verify binding retention"},
  {"name": "Asn_Deamidation_NS",    "pattern": r"N[ST]",
   "severity_cdr": "HIGH", "severity_fr": "LOW",
   "description": "Asn-Ser/Thr: moderate deamidation risk (t½ ~weeks)",
   "recommendation": "Monitor by forced deamidation study (pH 9, 40°C, 1 week)"},
  {"name": "Asn_Deamidation_other", "pattern": r"N[AHNQ]",
   "severity_cdr": "MEDIUM", "severity_fr": "LOW",
   "description": "Asn-X deamidation: lower risk motif",
   "recommendation": "Note for stability monitoring"},
  # Asp isomerization
  {"name": "Asp_Isomerization_DG",  "pattern": r"DG",
   "severity_cdr": "HIGH", "severity_fr": "MEDIUM",
   "description": "Asp-Gly: fastest Asp isomerization (succinimide formation)",
   "recommendation": "Mutate D→E (conservative) or avoid in CDR; check pH 4 stability"},
  {"name": "Asp_Isomerization_DS",  "pattern": r"D[ST]",
   "severity_cdr": "HIGH", "severity_fr": "LOW",
   "description": "Asp-Ser/Thr: moderate isomerization risk",
   "recommendation": "Assess by forced acidic stress (pH 4, 37°C)"},
  # Oxidation
  {"name": "Met_Oxidation",  "pattern": r"M",
   "severity_cdr": "HIGH", "severity_fr": "LOW",
   "description": "Methionine: oxidation risk (especially surface-exposed Met)",
   "recommendation": "Mutate M→L or M→V in CDRs; verify by forced oxidation (H₂O₂)"},
  {"name": "Trp_Oxidation",  "pattern": r"W",
   "severity_cdr": "HIGH", "severity_fr": "LOW",
   "description": "Tryptophan: photo-oxidation and stress oxidation risk",
   "recommendation": "Light-protect formulation; assess photostability (ICH Q1B)"},
  # N-linked glycosylation
  {"name": "N_Glycosylation_NxS_T", "pattern": r"N[^P][ST]",
   "severity_cdr": "HIGH", "severity_fr": "HIGH",
   "description": "N-X-S/T: N-linked glycosylation sequon (unintended variable-domain glycosylation)",
   "recommendation": "Mutate N→Q or S/T→A to remove; unintended glycosylation causes heterogeneity"},
  # Unpaired cysteine
  {"name": "Unpaired_Cys",   "pattern": r"C",
   "severity_cdr": "HIGH", "severity_fr": "MEDIUM",
   "description": "Free cysteine: risk of disulfide scrambling, aggregation, conjugation heterogeneity",
   "recommendation": "Verify paired in disulfide bond; mutate to Ser if unpaired"},
  # Proteolytic cleavage
  {"name": "Asp_Pro_Cleavage", "pattern": r"DP",
   "severity_cdr": "HIGH", "severity_fr": "MEDIUM",
   "description": "Asp-Pro: acid-labile peptide bond (risk at low pH purification steps)",
   "recommendation": "Avoid DP in CDRs; check stability at pH 3 (Protein A elution conditions)"},
  # RGD integrin-binding
  {"name": "RGD_Integrin",    "pattern": r"RGD",
   "severity_cdr": "HIGH", "severity_fr": "MEDIUM",
   "description": "RGD: integrin-binding motif → off-target binding risk",
   "recommendation": "Redesign if in CDR; RGD may cause polyreactivity"},
  # Lysine glycation
  {"name": "Lys_Glycation",   "pattern": r"K",
   "severity_cdr": "MEDIUM", "severity_fr": "LOW",
   "description": "Lysine: glycation risk (reductive amination with glucose in serum)",
   "recommendation": "Monitor CDR Lys by glycation MS; consider K→R substitution if problematic"},
]

def scan_liabilities(chain_data):
  """
  Scan a numbered antibody chain for all chemical liability motifs.
  CDR liabilities: HIGH severity. Framework liabilities: MEDIUM or LOW.
  """
  sequence = chain_data["sequence"]
  chain_obj = chain_data.get("chain_obj")
  chain_type = chain_data.get("chain_type", "unknown")

  # Build region map
  region_map = {}
  if chain_obj is not None:
    cdr_ranges = IMGT_CDR_POSITIONS.get(chain_type, {})
    pos_idx = 0
    for imgt_pos, aa in chain_obj:
      in_cdr = None
      for cdr_name, pos_range in cdr_ranges.items():
        if imgt_pos.number in pos_range:
          in_cdr = cdr_name
          break
      region_map[pos_idx] = in_cdr if in_cdr else "Framework"
      pos_idx += 1
  else:
    for i in range(len(sequence)):
      region_map[i] = "Unknown"

  liabilities = []
  for lib_def in LIABILITY_PATTERNS:
    pattern = lib_def["pattern"]
    for match in re.finditer(pattern, sequence):
      pos = match.start()
      region = region_map.get(pos, "Unknown")
      severity = lib_def["severity_cdr"] if "CDR" in str(region) else lib_def["severity_fr"]
      if severity == "LOW" and "CDR" not in str(region):
        continue
      liabilities.append(Liability(
        name=lib_def["name"],
        motif=match.group(),
        position=pos + 1,
        region=str(region),
        severity=severity,
        description=lib_def["description"],
        recommendation=lib_def["recommendation"],
      ))

  severity_order = {"HIGH": 0, "MEDIUM": 1, "LOW": 2}
  liabilities.sort(key=lambda x: (severity_order[x.severity], x.position))
  return liabilities
```

---

### Step 3: Physicochemical Profiling — TAP Five Metrics (Layer 2)

```python
import numpy as np

# Kyte-Doolittle hydrophobicity scale
KD_HYDROPHOBICITY = {
  "I": 4.5, "V": 4.2, "L": 3.8, "F": 2.8, "C": 2.5,
  "M": 1.9, "A": 1.8, "G": -0.4, "T": -0.7, "W": -0.9,
  "S": -0.8, "Y": -1.3, "P": -1.6, "H": -3.2, "E": -3.5,
  "Q": -3.5, "D": -3.5, "N": -3.5, "K": -3.9, "R": -4.5,
}
CHARGE_POS = {"R": 1.0, "K": 1.0, "H": 0.1}
CHARGE_NEG = {"D": -1.0, "E": -1.0}

# TAP CST reference distributions (Raybould et al. 2019, PNAS)
TAP_REFERENCE = {
  "CDR_total_length":  {"median": 44,  "p5": 34,  "p95": 55,  "unit": "aa"},
  "surface_hydrophob": {"median": 0.0,  "p5": -1.5, "p95": 2.5, "unit": "score"},
  "CDR_pos_charge":    {"median": 2.0,  "p5": 0.0,  "p95": 5.0, "unit": "charge"},
  "CDR_neg_charge":    {"median": -1.5, "p5": -4.0, "p95": 0.0, "unit": "charge"},
  "charge_asymmetry":  {"median": 0.0,  "p5": -2.0, "p95": 2.0, "unit": "ΔQ (H-L)"},
}

def estimate_pi(sequence):
  """Estimate pI using bisection charge-balance method."""
  pka = {"D": 3.65, "E": 4.25, "H": 6.00, "C": 8.18, "Y": 10.07, "K": 10.53, "R": 12.48}
  def net_charge_at_ph(ph):
    charge = 0.0
    for aa in sequence:
      if aa in pka:
        pk = pka[aa]
        charge -= 1/(1+10**(ph-pk)) if aa in "DE" else 1/(1+10**(pk-ph))
    charge += 1/(1+10**(9.0-ph))   # N-terminus
    charge -= 1/(1+10**(ph-2.0))   # C-terminus
    return charge
  lo, hi = 2.0, 14.0
  for _ in range(50):
    mid = (lo+hi)/2
    if net_charge_at_ph(mid) > 0: lo = mid
    else: hi = mid
  return round((lo+hi)/2, 2)

def compute_physicochemical_profile(chain_data, all_chains=None):
  """
  Compute the five TAP developability metrics (Raybould et al. PNAS 2019):
  1. CDR total length
  2. CDR surface hydrophobicity (KD sum)
  3. CDR positive charge
  4. CDR negative charge
  5. Charge asymmetry (VH - VL)
  """
  cdr_seq = chain_data.get("cdr_combined", chain_data["sequence"])
  cdr_length = chain_data.get("cdr_total_length", len(cdr_seq))
  hydrophob_score = sum(KD_HYDROPHOBICITY.get(aa, 0) for aa in cdr_seq)
  pos_charge = sum(CHARGE_POS.get(aa, 0) for aa in cdr_seq)
  neg_charge = sum(CHARGE_NEG.get(aa, 0) for aa in cdr_seq)

  charge_asym = 0.0
  if all_chains and len(all_chains) == 2:
    chain_list = list(all_chains.values())
    q1 = sum(CHARGE_POS.get(aa, 0) + CHARGE_NEG.get(aa, 0)
              for aa in chain_list[0].get("cdr_combined", ""))
    q2 = sum(CHARGE_POS.get(aa, 0) + CHARGE_NEG.get(aa, 0)
              for aa in chain_list[1].get("cdr_combined", ""))
    charge_asym = q1 - q2

  pi_estimate = estimate_pi(chain_data["sequence"])
  sequence = chain_data["sequence"]
  INSTAB_DIPEPTIDES = {"WW":1,"WC":1,"WM":1,"WH":1,"WD":1,"WF":1,"WK":1,
                       "FW":1,"CW":1,"MW":1,"AW":1,"DW":1,"HW":1}
  instab_count = sum(1 for i in range(len(sequence)-1) if sequence[i:i+2] in INSTAB_DIPEPTIDES)

  metrics = {
    "CDR_total_length":       cdr_length,
    "CDR_hydrophobicity":   round(hydrophob_score, 2),
    "CDR_positive_charge":  round(pos_charge, 2),
    "CDR_negative_charge":  round(neg_charge, 2),
    "charge_asymmetry_VH_VL": round(charge_asym, 2),
    "pI_estimate":           pi_estimate,
    "sequence_length":       len(sequence),
    "instability_dipeptide_count": instab_count,
  }

  tap_flags = {}
  tap_map = {
    "CDR_total_length":       "CDR_total_length",
    "CDR_hydrophobicity":    "surface_hydrophob",
    "CDR_positive_charge":   "CDR_pos_charge",
    "CDR_negative_charge":   "CDR_neg_charge",
    "charge_asymmetry_VH_VL": "charge_asymmetry",
  }
  for metric_key, ref_key in tap_map.items():
    ref = TAP_REFERENCE[ref_key]
    val = metrics[metric_key]
    within = ref["p5"] <= val <= ref["p95"]
    tap_flags[metric_key] = {
      "value": val, "within_CST": within,
      "flag": "OK" if within else "OUTSIDE CST",
      "reference": f"CST 90%: [{ref['p5']}, {ref['p95']}] {ref['unit']}",
    }

  print(f"\n--- TAP Metrics ---")
  for k, info in tap_flags.items():
    print(f"  {k}: {info['value']} — {info['flag']}")
  print(f"  pI: {pi_estimate}, length: {len(sequence)} aa")

  return {"metrics": metrics, "tap_flags": tap_flags}
```

---

### Step 4: Thera-SAbDab Benchmarking (Layer 3)

```python
import requests
import pandas as pd
from io import StringIO

THERASABDAB_URL = (
  "https://opig.stats.ox.ac.uk/webapps/sabdab-sabpred/static/downloads/"
  "TheraSAbDab_SeqStruc_OnlineDownload.csv"
)

def fetch_therasabdab(cache_path="therasabdab_cache.csv"):
  """Download all WHO-recognised therapeutic antibody sequences."""
  import os
  if os.path.exists(cache_path):
    return pd.read_csv(cache_path)
  print(f"Downloading Thera-SAbDab...")
  response = requests.get(THERASABDAB_URL, timeout=120)
  response.raise_for_status()
  df = pd.read_csv(StringIO(response.text))
  df.to_csv(cache_path, index=False)
  print(f"Downloaded {len(df)} therapeutic sequences")
  return df

def compute_sequence_identity(seq1, seq2):
  """Fast k-mer Jaccard similarity as proxy for sequence identity."""
  if not seq1 or not seq2: return 0.0
  if len(seq1) == len(seq2):
    return sum(a==b for a,b in zip(seq1,seq2)) / len(seq1)
  k = 3
  kmers1 = set(seq1[i:i+k] for i in range(len(seq1)-k+1))
  kmers2 = set(seq2[i:i+k] for i in range(len(seq2)-k+1))
  if not kmers1 or not kmers2: return 0.0
  return len(kmers1 & kmers2) / len(kmers1 | kmers2)

def benchmark_against_therapeutics(chain_data, therasabdab_df, top_n=5):
  """
  Find most similar approved/clinical-stage therapeutic antibodies.
  High similarity (>90%) → biosimilar territory
  Moderate (70-90%) → precedented structural space
  Low (<70%) → novel scaffold
  """
  if therasabdab_df is None:
    return {"status": "skipped", "reason": "download failed"}

  query_seq = chain_data["sequence"]
  vh_col = next((c for c in therasabdab_df.columns
                  if "heavy" in c.lower() or "vh" in c.lower()), None)
  if vh_col is None:
    return {"status": "skipped", "reason": "column unrecognised"}

  similarities = []
  for _, row in therasabdab_df.iterrows():
    ref_seq = str(row.get(vh_col, "")).strip()
    if len(ref_seq) < 50: continue
    sim = compute_sequence_identity(query_seq, ref_seq)
    similarities.append({
      "INN": row.get("INN", "Unknown"),
      "target": row.get("Target", "Unknown"),
      "clinical_stage": row.get("Clinical_Stage", "Unknown"),
      "sequence_identity": round(sim*100, 1),
    })

  similarities.sort(key=lambda x: x["sequence_identity"], reverse=True)
  top_matches = similarities[:top_n]
  best_sim = top_matches[0]["sequence_identity"]

  if best_sim > 90:   novelty = "LOW (potential IP/biosimilar issue)"
  elif best_sim > 70: novelty = "MODERATE (precedented structural space)"
  else:               novelty = "HIGH (novel sequence space)"

  print(f"\n--- Thera-SAbDab Benchmark ---")
  for m in top_matches:
    print(f"  {m['INN']} | {m['target']} | {m['clinical_stage']} | {m['sequence_identity']}%")
  print(f"  Novelty: {novelty}")

  return {
    "status": "success", "top_matches": top_matches,
    "best_identity": best_sim, "novelty": novelty,
  }
```

---

### Step 5: Developability Score and Traffic-Light Classification

```python
def compute_developability_score(liabilities, physchem, benchmark):
  """
  Composite score 0–100 (higher = lower risk).
  Weights: Chemical liabilities 40pts, TAP compliance 35pts, Benchmark 25pts.

  Traffic light:
    GREEN  score ≥ 75  → proceed, monitor
    AMBER  score 50–74 → engineer before advancing
    RED    score < 50  → significant risk, redesign recommended
  """
  # Liability penalty
  liability_penalty = 0
  for lib in liabilities:
    if "CDR" in str(lib.region):
      penalty = 8 if lib.severity == "HIGH" else 4
    else:
      penalty = 3 if lib.severity == "HIGH" else 1
    liability_penalty += penalty
  liability_score = max(0, 40 - liability_penalty)

  # TAP compliance
  tap_flags = physchem.get("tap_flags", {})
  n_within = sum(1 for f in tap_flags.values() if f.get("within_CST", False))
  tap_score = round(35 * n_within / (len(tap_flags) or 1))

  # Benchmark
  best_id = benchmark.get("best_identity", 0)
  if best_id > 95:   novelty_score = 10
  elif best_id > 80: novelty_score = 20
  elif best_id > 60: novelty_score = 25
  else:              novelty_score = 15

  total_score = min(100, max(0, liability_score + tap_score + novelty_score))

  if total_score >= 75:
    traffic_light = "GREEN (proceed)"
    recommendation = "Good developability profile. Advance to biophysical characterization."
  elif total_score >= 50:
    traffic_light = "AMBER (caution)"
    recommendation = "Moderate risk. Address flagged liabilities before lead selection."
  else:
    traffic_light = "RED (risk)"
    recommendation = "High developability risk. Engineering campaign required before advancing."

  print(f"\n{'='*55}")
  print(f" DEVELOPABILITY SCORE: {total_score}/100  [{traffic_light}]")
  print(f"   Liability: {liability_score}/40 | TAP: {tap_score}/35 | Benchmark: {novelty_score}/25")
  print(f"   {recommendation}")
  print(f"{'='*55}")

  return {
    "total_score": total_score, "traffic_light": traffic_light,
    "recommendation": recommendation,
    "breakdown": {"liability": liability_score, "tap": tap_score, "novelty": novelty_score},
  }
```

---

### Step 6: Visualization

```python
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import numpy as np

def plot_scorecard(scorecard, liabilities, physchem, chain_name, out_path="scorecard.png"):
  """
  Generate a three-panel scorecard PNG:
  1. Composite score gauge (semicircle with needle)
  2. TAP metrics bar chart (green=OK, red=WARN)
  3. Liability summary by region and severity
  """
  fig = plt.figure(figsize=(16, 10))
  fig.suptitle(f"Antibody Developability Scorecard — {chain_name}",
               fontsize=16, fontweight="bold", y=0.98)

  # Panel 1: Score Gauge
  ax1 = fig.add_subplot(1, 3, 1)
  score = scorecard["total_score"]
  color = "#2ecc71" if score >= 75 else "#f39c12" if score >= 50 else "#e74c3c"

  theta = np.linspace(np.pi, 0, 300)
  ax1.plot(np.cos(theta), np.sin(theta), "k-", lw=2)
  for start, end, col in [(0,50,"#e74c3c"), (50,75,"#f39c12"), (75,100,"#2ecc71")]:
    t = np.linspace(np.pi*(1-start/100), np.pi*(1-end/100), 50)
    ax1.fill_between(np.cos(t), 0, np.sin(t), alpha=0.3, color=col)

  angle = np.pi * (1 - score/100)
  ax1.annotate("", xy=(0.85*np.cos(angle), 0.85*np.sin(angle)),
               xytext=(0,0), arrowprops=dict(arrowstyle="-|>", color=color, lw=2))
  ax1.text(0, -0.15, f"{score}/100", ha="center", va="center", fontsize=28, fontweight="bold", color=color)
  tl = scorecard["traffic_light"]
  ax1.text(0, -0.32, tl, ha="center", fontsize=14, fontweight="bold", color=color)
  ax1.set_xlim(-1.2, 1.2); ax1.set_ylim(-0.5, 1.2); ax1.axis("off")
  ax1.set_title("Composite Score", fontsize=12, fontweight="bold")

  # Panel 2: TAP Metrics
  ax2 = fig.add_subplot(1, 3, 2)
  tap_flags = physchem.get("tap_flags", {})
  labels = {"CDR_total_length":"CDR Length","CDR_hydrophobicity":"Hydrophobicity",
            "CDR_positive_charge":"Pos. Charge","CDR_negative_charge":"Neg. Charge",
            "charge_asymmetry_VH_VL":"Charge Asym."}
  values = [v["value"] for v in tap_flags.values()]
  colors  = ["#2ecc71" if v["within_CST"] else "#e74c3c" for v in tap_flags.values()]
  y_pos   = range(len(tap_flags))
  ax2.barh(list(y_pos), values, color=colors, height=0.6, edgecolor="white")
  ax2.set_yticks(list(y_pos))
  ax2.set_yticklabels([labels.get(k,k) for k in tap_flags.keys()], fontsize=10)
  ax2.axvline(0, color="k", lw=0.8)
  for i, (val, v) in enumerate(zip(values, tap_flags.values())):
    flag = "OK" if v["within_CST"] else "WARN"
    offset = val + 0.1 if val >= 0 else val - 0.1
    ax2.text(offset, i, f"[{flag}] {val:.1f}", va="center", fontsize=9)
  ax2.set_title("TAP Metrics (vs CST Reference)", fontsize=11, fontweight="bold")
  ax2.legend(handles=[mpatches.Patch(color="#2ecc71"), mpatches.Patch(color="#e74c3c")],
             labels=["Within CST 90%","Outside CST 90%"], fontsize=8)

  # Panel 3: Liability Summary
  ax3 = fig.add_subplot(1, 3, 3)
  if liabilities:
    from collections import defaultdict
    by_region = defaultdict(lambda: {"HIGH":0,"MEDIUM":0})
    for lib in liabilities:
      by_region[str(lib.region)][lib.severity] += 1
    regions = list(by_region.keys())
    x = np.arange(len(regions)); w = 0.35
    ax3.bar(x-w/2, [by_region[r]["HIGH"] for r in regions], w, label="HIGH", color="#e74c3c", alpha=0.85)
    ax3.bar(x+w/2, [by_region[r]["MEDIUM"] for r in regions], w, label="MEDIUM", color="#f39c12", alpha=0.85)
    ax3.set_xticks(x); ax3.set_xticklabels(regions, rotation=30, ha="right", fontsize=9)
    ax3.set_ylabel("Liability Count"); ax3.set_title("Liabilities by Region & Severity", fontsize=11, fontweight="bold")
    ax3.legend(fontsize=9)
  else:
    ax3.text(0.5, 0.5, "No significant\nliabilities detected",
             ha="center", va="center", fontsize=14, color="#2ecc71", transform=ax3.transAxes)
    ax3.axis("off")

  plt.tight_layout(rect=[0, 0, 1, 0.95])
  plt.savefig(out_path, dpi=150, bbox_inches="tight")
  plt.close()
  print(f"Scorecard saved: {out_path}")
```

---

### Step 7: Main Orchestration

```python
import json, os

def run_abdev(sequence_input, name="query", out_dir="abdev_results", skip_benchmark=False):
  """
  Full AbDev pipeline entry point.
  Returns: Complete assessment dict with all metrics, liabilities, and scorecard.
  """
  os.makedirs(out_dir, exist_ok=True)
  print(f"\n{'='*55}")
  print(f" AbDev — Antibody Developability Assessment  |  Candidate: {name}")
  print(f"{'='*55}")

  chains = parse_antibody_input(sequence_input)
  all_results = {}

  for chain_name, chain_data in chains.items():
    print(f"\n{'─'*55}")
    print(f" Processing: {chain_name} ({chain_data['chain_type']})")

    liabilities   = scan_liabilities(chain_data)
    physchem      = compute_physicochemical_profile(chain_data, all_chains=chains)

    if not skip_benchmark:
      therasabdab_df = fetch_therasabdab(cache_path=os.path.join(out_dir, "therasabdab_cache.csv"))
      benchmark      = benchmark_against_therapeutics(chain_data, therasabdab_df)
    else:
      benchmark = {"status": "skipped", "best_identity": 70, "novelty": "N/A"}

    scorecard = compute_developability_score(liabilities, physchem, benchmark)

    scorecard_path = os.path.join(out_dir, f"scorecard_{chain_name}.png")
    plot_scorecard(scorecard, liabilities, physchem,
                   chain_name=f"{name} — {chain_name}", out_path=scorecard_path)

    all_results[chain_name] = {
      "chain_type": chain_data["chain_type"],
      "sequence_length": chain_data["length"],
      "cdr_lengths": {k: len(v) for k, v in chain_data["cdrs"].items()},
      "liabilities": [
        {"name": l.name, "position": l.position, "region": l.region,
         "severity": l.severity, "recommendation": l.recommendation}
        for l in liabilities
      ],
      "physicochemical": physchem["metrics"],
      "tap_flags": {k: v["flag"] for k, v in physchem["tap_flags"].items()},
      "benchmark": benchmark,
      "scorecard": scorecard,
    }

  # Save JSON
  out_json = os.path.join(out_dir, f"{name}_results.json")
  with open(out_json, "w") as f:
    json.dump({k: {x: v for x, v in r.items() if x != "chain_obj"}
               for k, r in all_results.items()}, f, indent=2, default=str)
  print(f"\nResults saved: {out_json}")
  return all_results

# ─── DEMO ─────────────────────────────────────────────────────────────────
if __name__ == "__main__":
  TRASTUZUMAB_VH = (
    "EVQLVESGGGLVQPGGSLRLSCAASGFNIKDTYIHWVRQAPGKGLEWVARIYPTNGYTRYADSVKGRFTISADTSKNTAYLQM"
    "NSLRAEDTAVYYCSRWGGDGFYAMDYWGQGTLVTVSS"
  )
  TRASTUZUMAB_VL = (
    "DIQMTQSPSSLSASVGDRVTITCRASQDVNTAVAWYQQKPGKAPKLLIYSASFLYSGVPSRFSGSRSGTDFTLTISSLQPED"
    "FATYYCQQHYTTPPTFGQGTKVEIK"
  )
  results = run_abdev(
    sequence_input=f"{TRASTUZUMAB_VH}:{TRASTUZUMAB_VL}",
    name="Trastuzumab_demo",
    out_dir="abdev_results_trastuzumab",
    skip_benchmark=False,
  )
```

---

## Developability Metrics Reference

| Metric | HIGH Risk Threshold | Source |
|--------|--------------------|--------|
| CDR total length | > 55 aa | Raybould et al. 2019 |
| CDR hydrophobicity | > 2.5 (KD sum) | Raybould et al. 2019 |
| CDR positive charge | > 5.0 | Raybould et al. 2019 |
| CDR negative charge | < -4.0 | Raybould et al. 2019 |
| Charge asymmetry VH-VL | > ±2.0 | Raybould et al. 2019 |
| Asn-Gly (NG) in CDR | Any occurrence | Lu et al. mAbs 2019 |
| N-X-S/T glycosylation | Any in variable domain | Industry standard |
| Asp-Gly (DG) in CDR | Any occurrence | Industry standard |

---

## Adaptation

- **Nanobody / VHH:** Pass single VHH sequence; ANARCI auto-detects camelid heavy-chain type
- **scFv:** Pass VH:VL as single string with colon separator
- **Batch screening:** Loop over multiple sequences and collect scorecards for ranked triage
- **Custom liability panel:** Add entries to `LIABILITY_PATTERNS` for molecule-specific concerns
- **Offline mode:** Use `--skip-benchmark` to skip Thera-SAbDab download

---

## Dependencies

```
# Required Python packages
abnumber>=0.4.4       # ANARCI wrapper for antibody numbering (includes HMMER dep)
pandas>=1.5
numpy>=1.24
matplotlib>=3.7
requests>=2.28
scipy>=1.10
biopython>=1.81

# System (HMMER — required by abnumber/ANARCI)
# Ubuntu/Debian:
sudo apt-get install hmmer
# Or via conda:
conda install -c bioconda hmmer -y
```

Python 3.9+. CPU only. Typical runtime: < 60 seconds per sequence.

---

## References

1. Raybould, M.I.J. et al. (2019). Five computational developability guidelines for therapeutic antibody profiling. *PNAS*.
2. Lu, X. et al. (2019). Deamidation and isomerization liability analysis of 131 clinical-stage antibodies. *mAbs*.
3. Raybould, M.I.J. et al. (2020). Thera-SAbDab: the Therapeutic Structural Antibody Database. *Nucleic Acids Research*.
4. Dunbar, J. & Deane, C.M. (2016). ANARCI: antigen receptor numbering and receptor classification. *Bioinformatics*.
5. Sharma, V.K. et al. (2023). Blueprint for antibody biologics developability. *mAbs*.
