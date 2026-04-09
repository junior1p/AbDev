#!/usr/bin/env python3
"""
AbDev — Antibody Developability Assessment Pipeline
===================================================
Fully automated in-silico developability profiling for therapeutic antibodies
and nanobodies. Produces an industry-standard scorecard from a single
amino acid sequence.

Three assessment layers:
  Layer 1 — Chemical liability scanning (deamidation, oxidation, isomerization, ...)
  Layer 2 — TAP physicochemical profiling (Raybould et al. PNAS 2019)
  Layer 3 — Thera-SAbDab benchmarking against approved therapeutics

No GPU, no external API keys. Python 3.9+ only.

Usage:
  python abdev_pipeline.py                    # run demo (Trastuzumab)
  python abdev_pipeline.py --seq "EVQLV..."   # single sequence
  python abdev_pipeline.py --fasta seqs.fasta # batch mode

Output:
  abdev_results/<name>/
    <name>_scorecard.png   — visual one-page scorecard
    <name>_results.json    — machine-readable full assessment
    therasabdab_cache.csv  — cached approved therapeutic sequences
"""

import argparse
import json
import logging
import os
import re
import sys
import warnings
from dataclasses import dataclass, field
from typing import Optional

# ── Third-party imports ──────────────────────────────────────────────────────
import requests
import pandas as pd
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from io import StringIO

try:
    from abnumber import Chain
    _HAS_ABNUMBER = True
except ImportError:
    _HAS_ABNUMBER = False
    warnings.warn("abnumber not installed. ANARCI numbering disabled. "
                  "Install with: pip install abnumber anarci")

# ── Logging ─────────────────────────────────────────────────────────────────
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s - AbDev - %(levelname)s - %(message)s",
    datefmt="%Y-%m-%d %H:%M:%S",
)
log = logging.getLogger("AbDev")

# ═══════════════════════════════════════════════════════════════════════════════
# CONSTANTS & REFERENCE DATA
# ═══════════════════════════════════════════════════════════════════════════════

# IMGT CDR position ranges
IMGT_CDR_POSITIONS = {
    "H": {"CDR1": range(27, 39), "CDR2": range(56, 66), "CDR3": range(105, 118)},
    "L": {"CDR1": range(27, 39), "CDR2": range(56, 66), "CDR3": range(105, 118)},
    "K": {"CDR1": range(27, 39), "CDR2": range(56, 66), "CDR3": range(105, 118)},
}

# Kyte-Doolittle hydrophobicity scale
KD_HYDROPHOBICITY = {
    "I": 4.5, "V": 4.2, "L": 3.8, "F": 2.8, "C": 2.5,
    "M": 1.9, "A": 1.8, "G": -0.4, "T": -0.7, "W": -0.9,
    "S": -0.8, "Y": -1.3, "P": -1.6, "H": -3.2, "E": -3.5,
    "Q": -3.5, "D": -3.5, "N": -3.5, "K": -3.9, "R": -4.5,
}

# Charge at pH 7.4
CHARGE_POS = {"R": 1.0, "K": 1.0, "H": 0.1}
CHARGE_NEG = {"D": -1.0, "E": -1.0}

# pKa values for pI estimation
PKA = {"D": 3.65, "E": 4.25, "H": 6.00, "C": 8.18, "Y": 10.07, "K": 10.53, "R": 12.48}

# TAP CST reference distributions (Raybould et al. 2019, PNAS)
# [median, 5th percentile, 95th percentile] from 242 post-phase-I therapeutics
TAP_REFERENCE = {
    "CDR_total_length":   {"median": 44,  "p5": 34,  "p95": 55,  "unit": "aa"},
    "surface_hydrophob":  {"median": 0.0,  "p5": -1.5, "p95": 2.5, "unit": "score"},
    "CDR_pos_charge":     {"median": 2.0,  "p5": 0.0,  "p95": 5.0, "unit": "charge"},
    "CDR_neg_charge":     {"median": -1.5, "p5": -4.0, "p95": 0.0, "unit": "charge"},
    "charge_asymmetry":   {"median": 0.0,  "p5": -2.0, "p95": 2.0, "unit": "ΔQ (H-L)"},
}

# Instability dipeptides (Guruprasad 1990)
INSTAB_DIPEPTIDES = {
    "WW": 1, "WC": 1, "WM": 1, "WH": 1, "WD": 1, "WF": 1, "WK": 1,
    "FW": 1, "CW": 1, "MW": 1, "AW": 1, "DW": 1, "HW": 1,
}

# ── Liability motif definitions ───────────────────────────────────────────────
LIABILITY_PATTERNS = [
    # Deamidation
    {
        "name": "Asn_Deamidation_NG",
        "pattern": r"NG",
        "severity_cdr": "HIGH",
        "severity_fr": "MEDIUM",
        "description": "Asn-Gly: fastest deamidation motif (t½ ~1 day at pH 7.4)",
        "recommendation": "Mutate N→Q or N→A if in CDR paratope; verify binding retention",
    },
    {
        "name": "Asn_Deamidation_NS",
        "pattern": r"N[ST]",
        "severity_cdr": "HIGH",
        "severity_fr": "LOW",
        "description": "Asn-Ser/Thr: moderate deamidation risk (t½ ~weeks)",
        "recommendation": "Monitor by forced deamidation study (pH 9, 40°C, 1 week)",
    },
    {
        "name": "Asn_Deamidation_other",
        "pattern": r"N[AHNQ]",
        "severity_cdr": "MEDIUM",
        "severity_fr": "LOW",
        "description": "Asn-X deamidation: lower risk motif",
        "recommendation": "Note for stability monitoring",
    },
    # Isomerization
    {
        "name": "Asp_Isomerization_DG",
        "pattern": r"DG",
        "severity_cdr": "HIGH",
        "severity_fr": "MEDIUM",
        "description": "Asp-Gly: fastest Asp isomerization (succinimide formation)",
        "recommendation": "Mutate D→E (conservative) or avoid in CDR; check pH 4 stability",
    },
    {
        "name": "Asp_Isomerization_DS",
        "pattern": r"D[ST]",
        "severity_cdr": "HIGH",
        "severity_fr": "LOW",
        "description": "Asp-Ser/Thr: moderate isomerization risk",
        "recommendation": "Assess by forced acidic stress (pH 4, 37°C)",
    },
    # Oxidation
    {
        "name": "Met_Oxidation",
        "pattern": r"M",
        "severity_cdr": "HIGH",
        "severity_fr": "LOW",
        "description": "Methionine: oxidation risk (especially surface-exposed Met)",
        "recommendation": "Mutate M→L or M→V in CDRs; verify by forced oxidation (H₂O₂)",
    },
    {
        "name": "Trp_Oxidation",
        "pattern": r"W",
        "severity_cdr": "HIGH",
        "severity_fr": "LOW",
        "description": "Tryptophan: photo-oxidation and stress oxidation risk",
        "recommendation": "Light-protect formulation; assess photostability (ICH Q1B)",
    },
    # Glycosylation
    {
        "name": "N_Glycosylation_NxS_T",
        "pattern": r"N[^P][ST]",
        "severity_cdr": "HIGH",
        "severity_fr": "HIGH",
        "description": "N-X-S/T: N-linked glycosylation sequon (unintended variable-domain glycosylation)",
        "recommendation": "Mutate N→Q or S/T→A to remove; unintended glycosylation causes heterogeneity",
    },
    # Unpaired cysteine
    {
        "name": "Unpaired_Cys",
        "pattern": r"C",
        "severity_cdr": "HIGH",
        "severity_fr": "MEDIUM",
        "description": "Free cysteine: risk of disulfide scrambling, aggregation, conjugation heterogeneity",
        "recommendation": "Verify paired in disulfide bond; mutate to Ser if unpaired",
    },
    # Proteolytic cleavage
    {
        "name": "Asp_Pro_Cleavage",
        "pattern": r"DP",
        "severity_cdr": "HIGH",
        "severity_fr": "MEDIUM",
        "description": "Asp-Pro: acid-labile peptide bond (risk at low pH purification steps)",
        "recommendation": "Avoid DP in CDRs; check stability at pH 3 (Protein A elution conditions)",
    },
    # RGD integrin-binding
    {
        "name": "RGD_Integrin",
        "pattern": r"RGD",
        "severity_cdr": "HIGH",
        "severity_fr": "MEDIUM",
        "description": "RGD: integrin-binding motif → off-target binding risk",
        "recommendation": "Redesign if in CDR; RGD may cause polyreactivity",
    },
    # Lysine glycation
    {
        "name": "Lys_Glycation",
        "pattern": r"K",
        "severity_cdr": "MEDIUM",
        "severity_fr": "LOW",
        "description": "Lysine: glycation risk (reductive amination with glucose in serum)",
        "recommendation": "Monitor CDR Lys by glycation MS; consider K→R substitution if problematic",
    },
]


# ═══════════════════════════════════════════════════════════════════════════════
# DATA CLASSES
# ═══════════════════════════════════════════════════════════════════════════════

@dataclass
class Liability:
    name: str
    motif: str
    position: int       # 1-indexed
    region: str         # "CDR1", "CDR2", "CDR3", "Framework", "Unknown"
    severity: str       # "HIGH", "MEDIUM", "LOW"
    description: str
    recommendation: str


# ═══════════════════════════════════════════════════════════════════════════════
# STEP 1 — INPUT PARSING & ANTIBODY NUMBERING
# ═══════════════════════════════════════════════════════════════════════════════

def parse_antibody_input(input_str: str) -> dict:
    """
    Parse antibody sequence input. Accepts:
      - Single sequence string (VH or VHH/nanobody)
      - Two sequences separated by ':' (VH:VL)
      - FASTA format (single or two entries)
      - dict: {"VH": "...", "VL": "..."}

    Returns dict of chain annotations with IMGT numbering via abnumber.
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
            log.warning(f"{name} sequence too short ({len(seq)} aa) — skipping ANARCI")
            numbered_chains[name] = {
                "sequence": seq, "chain_type": "unknown",
                "chain_obj": None, "cdrs": {}, "cdr_combined": seq,
                "framework": seq, "length": len(seq), "cdr_total_length": 0,
            }
            continue

        if _HAS_ABNUMBER:
            try:
                chain = Chain(seq, scheme="imgt")
                chain_type = chain.chain_type  # "H", "L", "K", or "VHH"

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
                log.info(f"{name}: {ct_label}, {len(seq)} aa, "
                         f"CDR lengths: H1={len(cdrs.get('CDR1',''))}, "
                         f"H2={len(cdrs.get('CDR2',''))}, "
                         f"H3={len(cdrs.get('CDR3',''))}")

            except Exception as e:
                log.warning(f"ANARCI numbering failed for {name}: {e} — falling back to raw analysis")
                numbered_chains[name] = _fallback_analysis(name, seq)
        else:
            numbered_chains[name] = _fallback_analysis(name, seq)

    return numbered_chains


def _fallback_analysis(name: str, seq: str) -> dict:
    """Fallback when abnumber is not available."""
    log.warning("Falling back to raw sequence analysis (no CDR annotation)")
    return {
        "sequence": seq,
        "chain_type": "unknown",
        "chain_obj": None,
        "cdrs": {},
        "cdr_combined": seq,
        "framework": seq,
        "length": len(seq),
        "cdr_total_length": len(seq),
    }


# ═══════════════════════════════════════════════════════════════════════════════
# STEP 2 — CHEMICAL LIABILITY SCANNING
# ═══════════════════════════════════════════════════════════════════════════════

def scan_liabilities(chain_data: dict) -> list:
    """
    Scan a numbered antibody chain for all chemical liability motifs.
    CDR liabilities: HIGH severity (directly affect binding/activity)
    Framework liabilities: MEDIUM or LOW severity
    """
    sequence = chain_data["sequence"]
    cdrs = chain_data["cdrs"]
    chain_type = chain_data["chain_type"]
    chain_obj = chain_data.get("chain_obj")

    # Build region map
    region_map = {}
    if chain_obj is not None:
        cdr_ranges = IMGT_CDR_POSITIONS.get(chain_type, {})
        for imgt_pos, aa in chain_obj:
            in_cdr = None
            for cdr_name, pos_range in cdr_ranges.items():
                if imgt_pos.number in pos_range:
                    in_cdr = cdr_name
                    break
            # We track position by iterating the chain object's positions
        # Rebuild region_map by position index
        region_map = {}
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

            if "CDR" in str(region):
                severity = lib_def["severity_cdr"]
            else:
                severity = lib_def["severity_fr"]

            # Skip LOW severity in Framework
            if severity == "LOW" and "CDR" not in str(region):
                continue

            liabilities.append(Liability(
                name=lib_def["name"],
                motif=match.group(),
                position=pos + 1,  # 1-indexed
                region=str(region),
                severity=severity,
                description=lib_def["description"],
                recommendation=lib_def["recommendation"],
            ))

    # Sort: HIGH > MEDIUM > LOW, then by position
    severity_order = {"HIGH": 0, "MEDIUM": 1, "LOW": 2}
    liabilities.sort(key=lambda x: (severity_order[x.severity], x.position))
    return liabilities


def summarize_liabilities(liabilities: list) -> dict:
    """Count and classify liabilities by severity."""
    counts = {"HIGH": 0, "MEDIUM": 0, "LOW": 0}
    cdr_liabilities = []

    for lib in liabilities:
        counts[lib.severity] += 1
        if "CDR" in str(lib.region):
            cdr_liabilities.append(lib)

    log.info(f"--- Chemical Liability Summary ---")
    log.info(f"  HIGH severity:   {counts['HIGH']} (immediate engineering priority)")
    log.info(f"  MEDIUM severity: {counts['MEDIUM']}")
    log.info(f"  LOW severity:   {counts['LOW']}")
    log.info(f"  CDR liabilities: {len(cdr_liabilities)} (highest impact on potency)")

    if cdr_liabilities:
        log.info(f"\n  🔴 CDR Liabilities (must address):")
        for lib in cdr_liabilities[:8]:
            log.info(f"    [{lib.severity}] Pos {lib.position} {lib.region}: "
                    f"{lib.name} (motif: {lib.motif})")

    return {"counts": counts, "cdr_liabilities": cdr_liabilities, "total": len(liabilities)}


# ═══════════════════════════════════════════════════════════════════════════════
# STEP 3 — PHYSICOCHEMICAL PROFILING (TAP METRICS)
# ═══════════════════════════════════════════════════════════════════════════════

def net_charge_at_ph(sequence: str, ph: float) -> float:
    """Compute net charge of a sequence at a given pH."""
    charge = 0.0
    for aa in sequence:
        if aa in PKA:
            pka_val = PKA[aa]
            if aa in "DE":
                charge -= 1 / (1 + 10 ** (ph - pka_val))
            else:
                charge += 1 / (1 + 10 ** (pka_val - ph))
    # N-terminus (pKa ~9) and C-terminus (pKa ~2)
    charge += 1 / (1 + 10 ** (9.0 - ph))
    charge -= 1 / (1 + 10 ** (ph - 2.0))
    return charge


def estimate_pi(sequence: str) -> float:
    """Estimate isoelectric point using bisection."""
    lo, hi = 2.0, 14.0
    for _ in range(50):
        mid = (lo + hi) / 2
        if net_charge_at_ph(sequence, mid) > 0:
            lo = mid
        else:
            hi = mid
    return round((lo + hi) / 2, 2)


def compute_physicochemical_profile(chain_data: dict, all_chains: dict = None) -> dict:
    """
    Compute the five TAP developability metrics.
    Raybould et al. PNAS 2019 — five computational developability guidelines.
    """
    cdr_seq = chain_data.get("cdr_combined", chain_data["sequence"])
    cdrs = chain_data.get("cdrs", {})

    # 1. CDR total length
    cdr_length = chain_data.get("cdr_total_length", len(cdr_seq))

    # 2. CDR hydrophobicity (sum of KD scores for CDR residues)
    hydrophob_score = sum(KD_HYDROPHOBICITY.get(aa, 0) for aa in cdr_seq)

    # 3. CDR positive charge
    pos_charge = sum(CHARGE_POS.get(aa, 0) for aa in cdr_seq)

    # 4. CDR negative charge
    neg_charge = sum(CHARGE_NEG.get(aa, 0) for aa in cdr_seq)

    # 5. Charge asymmetry (VH - VL)
    charge_asym = 0.0
    if all_chains and len(all_chains) == 2:
        chain_list = list(all_chains.values())
        q1 = sum(CHARGE_POS.get(aa, 0) + CHARGE_NEG.get(aa, 0)
                 for aa in chain_list[0].get("cdr_combined", ""))
        q2 = sum(CHARGE_POS.get(aa, 0) + CHARGE_NEG.get(aa, 0)
                 for aa in chain_list[1].get("cdr_combined", ""))
        charge_asym = q1 - q2

    # pI estimate
    pi_estimate = estimate_pi(chain_data["sequence"])

    # Instability dipeptide count
    sequence = chain_data["sequence"]
    instab_count = sum(
        1 for i in range(len(sequence) - 1)
        if sequence[i:i+2] in INSTAB_DIPEPTIDES
    )

    metrics = {
        "CDR_total_length":         cdr_length,
        "CDR_hydrophobicity":       round(hydrophob_score, 2),
        "CDR_positive_charge":      round(pos_charge, 2),
        "CDR_negative_charge":     round(neg_charge, 2),
        "charge_asymmetry_VH_VL":   round(charge_asym, 2),
        "pI_estimate":              pi_estimate,
        "sequence_length":          len(sequence),
        "instability_dipeptide_count": instab_count,
    }

    # TAP flag: within CST 90% reference interval?
    tap_flags = {}
    tap_map = {
        "CDR_total_length":       "CDR_total_length",
        "CDR_hydrophobicity":     "surface_hydrophob",
        "CDR_positive_charge":    "CDR_pos_charge",
        "CDR_negative_charge":    "CDR_neg_charge",
        "charge_asymmetry_VH_VL": "charge_asymmetry",
    }

    for metric_key, ref_key in tap_map.items():
        ref = TAP_REFERENCE[ref_key]
        val = metrics[metric_key]
        within = ref["p5"] <= val <= ref["p95"]
        tap_flags[metric_key] = {
            "value":       val,
            "within_CST":  within,
            "flag":        "✅ OK" if within else "⚠️ OUTSIDE CST",
            "reference":   f"CST 90%: [{ref['p5']}, {ref['p95']}] {ref['unit']}",
        }

    log.info(f"\n--- Physicochemical Profile (TAP Metrics) ---")
    log.info(f"  {'Metric':<30} {'Value':>8}  {'Flag':<20}  {'Reference'}")
    log.info(f"  {'-'*75}")
    for k, info in tap_flags.items():
        log.info(f"  {k:<30} {info['value']:>8.2f}  {info['flag']:<20}  {info['reference']}")
    log.info(f"  pI (estimate): {pi_estimate}")
    log.info(f"  Sequence length: {len(sequence)} aa")

    return {"metrics": metrics, "tap_flags": tap_flags}


# ═══════════════════════════════════════════════════════════════════════════════
# STEP 4 — THERA-SAbDab BENCHMARKING
# ═══════════════════════════════════════════════════════════════════════════════

THERASABDAB_URL = "https://opig.stats.ox.ac.uk/webapps/sabdab-sabpred/static/downloads/TheraSAbDab_SeqStruc_OnlineDownload.csv"


def fetch_therasabdab(cache_path: str = "therasabdab_cache.csv") -> Optional[pd.DataFrame]:
    """Download all WHO-recognised therapeutic antibody sequences from Thera-SAbDab."""
    import os
    if os.path.exists(cache_path):
        log.info(f"Loading cached Thera-SAbDab sequences from {cache_path}")
        return pd.read_csv(cache_path)

    log.info(f"Downloading Thera-SAbDab from {THERASABDAB_URL}...")
    try:
        response = requests.get(THERASABDAB_URL, timeout=120)
        response.raise_for_status()
        df = pd.read_csv(StringIO(response.text))
        df.to_csv(cache_path, index=False)
        log.info(f"Downloaded {len(df)} therapeutic sequences to {cache_path}")
        return df
    except Exception as e:
        log.warning(f"Thera-SAbDab download failed: {e} — benchmarking skipped")
        return None


def compute_sequence_identity(seq1: str, seq2: str) -> float:
    """Fast pairwise sequence identity using k-mer Jaccard similarity."""
    if not seq1 or not seq2:
        return 0.0
    if len(seq1) == len(seq2):
        matches = sum(a == b for a, b in zip(seq1, seq2))
        return matches / len(seq1)
    k = 3
    kmers1 = set(seq1[i:i+k] for i in range(len(seq1) - k + 1))
    kmers2 = set(seq2[i:i+k] for i in range(len(seq2) - k + 1))
    if not kmers1 or not kmers2:
        return 0.0
    intersection = len(kmers1 & kmers2)
    union = len(kmers1 | kmers2)
    return intersection / union if union > 0 else 0.0


def benchmark_against_therapeutics(
    chain_data: dict,
    therasabdab_df: pd.DataFrame,
    top_n: int = 5,
) -> dict:
    """Find most similar approved/clinical-stage therapeutic antibodies."""
    if therasabdab_df is None:
        return {"status": "skipped", "reason": "Thera-SAbDab download failed"}

    query_seq = chain_data["sequence"]

    # Find the sequence column
    vh_col = next((
        c for c in therasabdab_df.columns
        if "heavy" in c.lower() or "vh" in c.lower() or "sequence" in c.lower()
    ), None)
    if vh_col is None:
        log.warning("Could not identify sequence column in Thera-SAbDab CSV")
        return {"status": "skipped", "reason": "Column format unrecognised"}

    # Compute similarity
    similarities = []
    for _, row in therasabdab_df.iterrows():
        ref_seq = str(row.get(vh_col, "")).strip()
        if len(ref_seq) < 50:
            continue
        sim = compute_sequence_identity(query_seq, ref_seq)
        inn = row.get("INN", row.get("Therapeutic", row.get("Name", "Unknown")))
        target = row.get("Target", row.get("Antigen", "Unknown"))
        stage = row.get("Clinical_Stage", row.get("Stage", "Unknown"))
        similarities.append({
            "INN":             inn,
            "target":          target,
            "clinical_stage":  stage,
            "sequence_identity": round(sim * 100, 1),
        })

    if not similarities:
        return {"status": "no_matches", "results": []}

    similarities.sort(key=lambda x: x["sequence_identity"], reverse=True)
    top_matches = similarities[:top_n]
    best_sim = top_matches[0]["sequence_identity"]

    if best_sim > 90:
        novelty = "LOW — very similar to approved therapeutic (potential IP/biosimilar issue)"
    elif best_sim > 70:
        novelty = "MODERATE — precedented structural space"
    else:
        novelty = "HIGH — novel sequence space"

    log.info(f"\n--- Thera-SAbDab Benchmarking (Top {top_n} Similar Therapeutics) ---")
    log.info(f"  {'INN':<25} {'Target':<20} {'Stage':<20} {'%Identity'}")
    log.info(f"  {'-'*75}")
    for m in top_matches:
        log.info(f"  {str(m['INN']):<25} {str(m['target']):<20} "
                 f"{str(m['clinical_stage']):<20} {m['sequence_identity']:.1f}%")
    log.info(f"\n  Novelty assessment: {novelty}")

    return {
        "status":          "success",
        "top_matches":     top_matches,
        "best_identity":   best_sim,
        "novelty":         novelty,
    }


# ═══════════════════════════════════════════════════════════════════════════════
# STEP 5 — DEVELOPABILITY SCORE & TRAFFIC LIGHT
# ═══════════════════════════════════════════════════════════════════════════════

def compute_developability_score(
    liabilities: list,
    physchem: dict,
    benchmark: dict,
) -> dict:
    """
    Compute composite developability score (0–100, higher = lower risk).

    Scoring weights:
      - Chemical liabilities:  40 pts  (most impactful on CMC)
      - TAP metrics compliance:  35 pts  (stability/manufacturability)
      - Novelty/benchmark:       25 pts  (IP/regulatory context)

    Traffic light:
      🟢 GREEN  score ≥ 75  → proceed, monitor
      🟡 AMBER  score 50–74  → engineer before advancing
      🔴 RED    score < 50   → significant risk, redesign recommended
    """
    # ── Liability penalty ──────────────────────────────────────────────────
    liability_penalty = 0
    for lib in liabilities:
        if "CDR" in str(lib.region):
            penalty = 8 if lib.severity == "HIGH" else 4
        else:
            penalty = 3 if lib.severity == "HIGH" else 1
        liability_penalty += penalty
    liability_score = max(0, 40 - liability_penalty)

    # ── TAP compliance score ────────────────────────────────────────────────
    tap_flags = physchem.get("tap_flags", {})
    n_within = sum(1 for f in tap_flags.values() if f.get("within_CST", False))
    n_total = len(tap_flags) or 1
    tap_score = round(35 * n_within / n_total)

    # ── Novelty/benchmark score ───────────────────────────────────────────
    best_id = benchmark.get("best_identity", 0)
    if best_id > 95:
        novelty_score = 10
    elif best_id > 80:
        novelty_score = 20
    elif best_id > 60:
        novelty_score = 25
    else:
        novelty_score = 15

    total_score = liability_score + tap_score + novelty_score
    total_score = min(100, max(0, total_score))

    if total_score >= 75:
        traffic_light = "🟢 GREEN"
        recommendation = "Good developability profile. Advance to biophysical characterization."
    elif total_score >= 50:
        traffic_light = "🟡 AMBER"
        recommendation = "Moderate risk. Address flagged liabilities before lead selection."
    else:
        traffic_light = "🔴 RED"
        recommendation = "High developability risk. Engineering campaign required before advancing."

    log.info(f"\n{'='*55}")
    log.info(f" DEVELOPABILITY SCORECARD")
    log.info(f"{'='*55}")
    log.info(f"  Liability score:      {liability_score:3d}/40")
    log.info(f"  TAP metrics score:   {tap_score:3d}/35")
    log.info(f"  Benchmark score:     {novelty_score:3d}/25")
    log.info(f"  {'─'*40}")
    log.info(f"  TOTAL SCORE: {total_score:3d}/100  {traffic_light}")
    log.info(f"\n  Recommendation: {recommendation}")
    log.info(f"{'='*55}")

    return {
        "total_score":    total_score,
        "traffic_light":  traffic_light,
        "recommendation": recommendation,
        "breakdown": {
            "liability": liability_score,
            "tap":        tap_score,
            "novelty":    novelty_score,
        },
    }


# ═══════════════════════════════════════════════════════════════════════════════
# STEP 6 — VISUALIZATION
# ═══════════════════════════════════════════════════════════════════════════════

def plot_scorecard(
    scorecard: dict,
    liabilities: list,
    physchem: dict,
    chain_name: str,
    out_path: str = "scorecard.png",
) -> str:
    """Generate a one-page visual scorecard (three-panel PNG)."""
    fig = plt.figure(figsize=(16, 10))
    fig.suptitle(f"Antibody Developability Scorecard — {chain_name}",
                 fontsize=16, fontweight="bold", y=0.98)

    # ── Panel 1: Score Gauge ────────────────────────────────────────────────
    ax1 = fig.add_subplot(1, 3, 1)
    score = scorecard["total_score"]
    color = "#2ecc71" if score >= 75 else "#f39c12" if score >= 50 else "#e74c3c"

    theta = np.linspace(np.pi, 0, 300)
    ax1.plot(np.cos(theta), np.sin(theta), "k-", linewidth=2)

    for start, end, col, label in [
        (0, 50, "#e74c3c", "RED (<50)"),
        (50, 75, "#f39c12", "AMBER (50-74)"),
        (75, 100, "#2ecc71", "GREEN (≥75)"),
    ]:
        t_start = np.pi * (1 - start / 100)
        t_end   = np.pi * (1 - end / 100)
        t = np.linspace(t_start, t_end, 50)
        ax1.fill_between(np.cos(t), 0, np.sin(t), alpha=0.3,
                         color=col.replace("80", ""))

    angle = np.pi * (1 - score / 100)
    ax1.annotate("", xy=(0.85 * np.cos(angle), 0.85 * np.sin(angle)),
                 xytext=(0, 0),
                 arrowprops=dict(arrowstyle="-|>", color=color, lw=2))
    ax1.text(0, -0.15, f"{score}/100", ha="center", va="center",
             fontsize=28, fontweight="bold", color=color)
    # Replace emoji traffic light with color-coded text
    tl = scorecard["traffic_light"]
    if "GREEN" in tl:
        tl_text = "GREEN (proceed)"
    elif "AMBER" in tl:
        tl_text = "AMBER (caution)"
    else:
        tl_text = "RED (risk)"
    ax1.text(0, -0.32, tl_text, ha="center",
             fontsize=14, fontweight="bold", color=color)
    ax1.set_xlim(-1.2, 1.2)
    ax1.set_ylim(-0.5, 1.2)
    ax1.axis("off")
    ax1.set_title("Composite Score", fontsize=12, fontweight="bold")

    # ── Panel 2: TAP Metrics Bars ──────────────────────────────────────────
    ax2 = fig.add_subplot(1, 3, 2)
    tap_flags = physchem.get("tap_flags", {})
    metric_labels = {
        "CDR_total_length":       "CDR Length",
        "CDR_hydrophobicity":    "Hydrophobicity",
        "CDR_positive_charge":    "Pos. Charge",
        "CDR_negative_charge":   "Neg. Charge",
        "charge_asymmetry_VH_VL": "Charge Asym.",
    }
    values = [v["value"] for v in tap_flags.values()]
    colors_bar = ["#2ecc71" if v["within_CST"] else "#e74c3c" for v in tap_flags.values()]
    status_labels = ["✅" if v["within_CST"] else "⚠️" for v in tap_flags.values()]

    y_pos = range(len(tap_flags))
    bars = ax2.barh(list(y_pos), values, color=colors_bar,
                    edgecolor="white", linewidth=0.5, height=0.6)
    ax2.set_yticks(list(y_pos))
    ax2.set_yticklabels([metric_labels.get(k, k) for k in tap_flags.keys()], fontsize=10)
    ax2.axvline(0, color="black", linewidth=0.8)
    for i, (val, status) in enumerate(zip(values, status_labels)):
        offset = val + 0.1 if val >= 0 else val - 0.1
        # Replace emoji with text to avoid DejaVu Sans glyph issues
        status_text = "OK" if status == "✅" else "WARN"
        ax2.text(offset, i, f"[{status_text}] {val:.1f}", va="center", fontsize=9)
    ax2.set_title("TAP Physicochemical Metrics\n(vs. CST Reference)", fontsize=11, fontweight="bold")
    ax2.set_xlabel("Score / Length", fontsize=9)
    green_p = mpatches.Patch(color="#2ecc71", label="Within CST 90%")
    red_p    = mpatches.Patch(color="#e74c3c", label="Outside CST 90%")
    ax2.legend(handles=[green_p, red_p], fontsize=8, loc="lower right")

    # ── Panel 3: Liability Summary ─────────────────────────────────────────
    ax3 = fig.add_subplot(1, 3, 3)
    if liabilities:
        from collections import defaultdict
        by_region = defaultdict(lambda: {"HIGH": 0, "MEDIUM": 0})
        for lib in liabilities:
            region = str(lib.region) if lib.region else "Unknown"
            by_region[region][lib.severity] += 1

        regions = list(by_region.keys())
        high_counts = [by_region[r]["HIGH"] for r in regions]
        med_counts  = [by_region[r]["MEDIUM"] for r in regions]
        x = np.arange(len(regions))
        w = 0.35
        ax3.bar(x - w/2, high_counts, w, label="HIGH",   color="#e74c3c", alpha=0.85)
        ax3.bar(x + w/2, med_counts,  w, label="MEDIUM",  color="#f39c12", alpha=0.85)
        ax3.set_xticks(x)
        ax3.set_xticklabels(regions, rotation=30, ha="right", fontsize=9)
        ax3.set_ylabel("Liability Count")
        ax3.set_title("Chemical Liabilities\nby Region & Severity", fontsize=11, fontweight="bold")
        ax3.legend(fontsize=9)
    else:
        ax3.text(0.5, 0.5, "No significant\nliabilities detected",
                 ha="center", va="center", fontsize=14, color="#2ecc71",
                 transform=ax3.transAxes)
        ax3.set_title("Chemical Liabilities", fontsize=11, fontweight="bold")
        ax3.axis("off")

    plt.tight_layout(rect=[0, 0, 1, 0.95])
    plt.savefig(out_path, dpi=150, bbox_inches="tight")
    plt.close()
    log.info(f"Scorecard saved: {out_path}")
    return out_path


# ═══════════════════════════════════════════════════════════════════════════════
# STEP 7 — MAIN ORCHESTRATION
# ═══════════════════════════════════════════════════════════════════════════════

def run_abdev(
    sequence_input,
    name: str = "query",
    out_dir: str = "abdev_results",
    skip_benchmark: bool = False,
) -> dict:
    """
    Full AbDev pipeline entry point.

    Args:
        sequence_input: Antibody sequence(s) — see parse_antibody_input()
        name: Label for the antibody candidate
        out_dir: Output directory
        skip_benchmark: Skip Thera-SAbDab download (offline mode)

    Returns:
        Complete assessment dict with all metrics, liabilities, and scorecard
    """
    os.makedirs(out_dir, exist_ok=True)

    log.info(f"\n{'='*55}")
    log.info(f" AbDev — Antibody Developability Assessment")
    log.info(f" Candidate: {name}")
    log.info(f"{'='*55}")

    # Step 1: Parse and number
    chains = parse_antibody_input(sequence_input)
    all_results = {}

    for chain_name, chain_data in chains.items():
        log.info(f"\n{'─'*55}")
        log.info(f" Processing: {chain_name} ({chain_data['chain_type']})")

        # Step 2: Liability scanning
        liabilities = scan_liabilities(chain_data)
        liability_summary = summarize_liabilities(liabilities)

        # Step 3: Physicochemical profiling
        physchem = compute_physicochemical_profile(chain_data, all_chains=chains)

        # Step 4: Thera-SAbDab benchmarking
        if not skip_benchmark:
            therasabdab_df = fetch_therasabdab(
                cache_path=os.path.join(out_dir, "therasabdab_cache.csv")
            )
            benchmark = benchmark_against_therapeutics(chain_data, therasabdab_df)
        else:
            benchmark = {"status": "skipped", "best_identity": 70, "novelty": "N/A"}

        # Step 5: Scorecard
        scorecard = compute_developability_score(liabilities, physchem, benchmark)

        # Step 6: Visualization
        scorecard_path = os.path.join(out_dir, f"scorecard_{chain_name}.png")
        plot_scorecard(
            scorecard, liabilities, physchem,
            chain_name=f"{name} — {chain_name}",
            out_path=scorecard_path,
        )

        all_results[chain_name] = {
            "chain_type":    chain_data["chain_type"],
            "sequence":      chain_data["sequence"],
            "sequence_length": chain_data["length"],
            "cdr_lengths":   {k: len(v) for k, v in chain_data["cdrs"].items()},
            "liabilities": [
                {
                    "name":          l.name,
                    "position":      l.position,
                    "region":        l.region,
                    "severity":      l.severity,
                    "description":   l.description,
                    "recommendation": l.recommendation,
                }
                for l in liabilities
            ],
            "liability_summary": liability_summary["counts"],
            "physicochemical":   physchem["metrics"],
            "tap_flags":         {k: v["flag"] for k, v in physchem["tap_flags"].items()},
            "benchmark":         benchmark,
            "scorecard":         scorecard,
        }

    # Save JSON
    safe_results = {}
    for k, v in all_results.items():
        safe_results[k] = {key: val for key, val in v.items() if key != "chain_obj"}

    out_json = os.path.join(out_dir, f"{name}_results.json")
    with open(out_json, "w") as f:
        json.dump(safe_results, f, indent=2, default=str)
    log.info(f"\nFull results saved: {out_json}")

    return all_results


# ═══════════════════════════════════════════════════════════════════════════════
# CLI & DEMO
# ═══════════════════════════════════════════════════════════════════════════════

if __name__ == "__main__":
    # Trastuzumab (Herceptin) — gold-standard FDA-approved anti-HER2 antibody
    TRASTUZUMAB_VH = (
        "EVQLVESGGGLVQPGGSLRLSCAASGFNIKDTYIHWVRQAPGKGLEWVARIYPTNGYTRYADSVKGRFTISADTSKNTAYLQM"
        "NSLRAEDTAVYYCSRWGGDGFYAMDYWGQGTLVTVSS"
    )
    TRASTUZUMAB_VL = (
        "DIQMTQSPSSLSASVGDRVTITCRASQDVNTAVAWYQQKPGKAPKLLIYSASFLYSGVPSRFSGSRSGTDFTLTISSLQPED"
        "FATYYCQQHYTTPPTFGQGTKVEIK"
    )

    parser = argparse.ArgumentParser(description="AbDev — Antibody Developability Assessment")
    parser.add_argument("--seq",    help="Single antibody sequence (VH or VHH)")
    parser.add_argument("--vhvl",   help="VH:VL pair separated by colon")
    parser.add_argument("--fasta",  help="FASTA file path")
    parser.add_argument("--name",   default="AbDev_query", help="Sample name")
    parser.add_argument("--out",    default="abdev_results", help="Output directory")
    parser.add_argument("--skip-benchmark", action="store_true", help="Skip Thera-SAbDab benchmarking")
    args = parser.parse_args()

    if args.fasta:
        with open(args.fasta) as f:
            sequence_input = f.read()
        name = os.path.splitext(os.path.basename(args.fasta))[0]
    elif args.vhvl:
        sequence_input = args.vhvl
        name = args.name
    elif args.seq:
        sequence_input = args.seq
        name = args.name
    else:
        # Demo: Trastuzumab
        sequence_input = f"{TRASTUZUMAB_VH}:{TRASTUZUMAB_VL}"
        name = "Trastuzumab_demo"
        args.out = "abdev_results_trastuzumab"
        log.info("Running demo with Trastuzumab (Herceptin) VH+VL sequences")

    results = run_abdev(
        sequence_input=sequence_input,
        name=name,
        out_dir=args.out,
        skip_benchmark=args.skip_benchmark,
    )

    # Print summary
    for chain_name, result in results.items():
        sc = result["scorecard"]
        print(f"\n{'='*55}")
        print(f" ✅ AbDev Complete: {name}")
        print(f"    Total Score: {sc['total_score']}/100  {sc['traffic_light']}")
        print(f"    Recommendation: {sc['recommendation']}")
        print(f"{'='*55}")
