#!/usr/bin/env python3
"""
AbDev Flask API Server
======================
REST API for the AbDev Antibody Developability Assessment Pipeline.

Endpoints:
  GET  /health          → health check
  POST /api/evaluate    → run full AbDev assessment
  GET  /                → web interface

Usage:
  python abdev_api.py                    # development (port 5001)
  gunicorn abdev_api:app -b 0.0.0.0:5001  # production
"""

import base64
import io
import logging
import os
import sys

from flask import (Flask, jsonify, request, render_template,
                   send_file, abort)
from werkzeug.exceptions import BadRequest

# Add parent dir to path so we can import abdev_pipeline
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from abdev_pipeline import run_abdev

# ── Logging ─────────────────────────────────────────────────────────────────
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s - AbDevAPI - %(levelname)s - %(message)s",
)
log = logging.getLogger("AbDevAPI")

# ── App ─────────────────────────────────────────────────────────────────────
app = Flask(__name__, template_folder="templates")
app.config["MAX_CONTENT_LENGTH"] = 16 * 1024 * 1024  # 16 MB

OUTPUT_DIR = os.environ.get("ABDEV_OUTPUT_DIR", "abdev_api_results")
os.makedirs(OUTPUT_DIR, exist_ok=True)


# ── Routes ──────────────────────────────────────────────────────────────────

@app.route("/")
def index():
    return render_template("index.html")


@app.route("/health")
def health():
    return jsonify({"status": "healthy", "service": "AbDev API", "version": "1.0.0"})


@app.route("/api/evaluate", methods=["POST"])
def api_evaluate():
    """
    Run full AbDev developability assessment.

    Request body (JSON):
    {
        "sequence": "EVQLVESGGGLVQPGGSLRLSCAAS...:DIQMTQSPSSLSASVGDRVTIT...",  // required
        "name":    "my_antibody",                                            // optional
        "skip_benchmark": false                                               // optional
    }

    sequence format:
      - Single VH/VHH:  "EVQLVESGGGLVQPGGSLRLSCAASGFNIKDTYIH..."
      - VH:VL pair:      "EVQLVESGGGLVQPGGSLRLSCAAS...:DIQMTQSPSSLSASVGDRVTIT..."
      - FASTA:           ">seq1\nEVQLVES...\n>seq2\nDIQMTQS..."

    Response (JSON):
    {
        "success": true,
        "name": "my_antibody",
        "results": {
            "Chain1": { "scorecard": { "total_score": 76, ... }, ... },
            "Chain2": { ... }
        },
        "scorecard_images": {
            "Chain1": "data:image/png;base64,...",
            "Chain2": "data:image/png;base64,..."
        }
    }
    """
    try:
        data = request.get_json()
    except Exception:
        return BadRequest("Invalid JSON body")

    if not data:
        return BadRequest("No JSON data provided")

    sequence = data.get("sequence", "").strip()
    if not sequence:
        return BadRequest("'sequence' field is required")

    name = data.get("name", "AbDev_query").strip()
    skip_benchmark = bool(data.get("skip_benchmark", False))

    log.info(f"Evaluating antibody: {name}")
    log.info(f"  Sequence length: {len(sequence)} chars")
    log.info(f"  skip_benchmark:  {skip_benchmark}")

    # Run pipeline in a dedicated output dir per request
    import time
    run_id = f"{name}_{int(time.time())}"
    run_out_dir = os.path.join(OUTPUT_DIR, run_id)

    try:
        results = run_abdev(
            sequence_input=sequence,
            name=name,
            out_dir=run_out_dir,
            skip_benchmark=skip_benchmark,
        )
    except Exception as e:
        log.error(f"Pipeline error: {e}")
        import traceback
        traceback.print_exc()
        return jsonify({
            "success": False,
            "error": str(e),
            "name": name,
        }), 500

    # Encode scorecard PNGs as base64
    scorecard_images = {}
    for chain_name in results:
        img_path = os.path.join(run_out_dir, f"scorecard_{chain_name}.png")
        if os.path.exists(img_path):
            with open(img_path, "rb") as f:
                img_b64 = base64.b64encode(f.read()).decode("utf-8")
            scorecard_images[chain_name] = f"data:image/png;base64,{img_b64}"
        else:
            scorecard_images[chain_name] = None

    # Build serialisable results dict (no chain_obj)
    serialised = {}
    for chain_name, result in results.items():
        serialised[chain_name] = {k: v for k, v in result.items() if k != "chain_obj"}

    log.info(f"Assessment complete: {name} — "
             f"scores: {[r['scorecard']['total_score'] for r in serialised.values()]}")

    return jsonify({
        "success": True,
        "name": name,
        "results": serialised,
        "scorecard_images": scorecard_images,
    })


@app.route("/api/demo", methods=["GET"])
def api_demo():
    """
    Run the built-in Trastuzumab demo (no request body needed).
    Useful for testing the API without a sequence.
    """
    TRASTUZUMAB_VH = (
        "EVQLVESGGGLVQPGGSLRLSCAASGFNIKDTYIHWVRQAPGKGLEWVARIYPTNGYTRYADSVKGRFTISADTSKNTAYLQM"
        "NSLRAEDTAVYYCSRWGGDGFYAMDYWGQGTLVTVSS"
    )
    TRASTUZUMAB_VL = (
        "DIQMTQSPSSLSASVGDRVTITCRASQDVNTAVAWYQQKPGKAPKLLIYSASFLYSGVPSRFSGSRSGTDFTLTISSLQPED"
        "FATYYCQQHYTTPPTFGQGTKVEIK"
    )
    sequence = f"{TRASTUZUMAB_VH}:{TRASTUZUMAB_VL}"
    name = "Trastuzumab_demo"

    import time
    run_id = f"demo_{int(time.time())}"
    run_out_dir = os.path.join(OUTPUT_DIR, run_id)

    results = run_abdev(
        sequence_input=sequence,
        name=name,
        out_dir=run_out_dir,
        skip_benchmark=False,
    )

    scorecard_images = {}
    for chain_name in results:
        img_path = os.path.join(run_out_dir, f"scorecard_{chain_name}.png")
        if os.path.exists(img_path):
            with open(img_path, "rb") as f:
                img_b64 = base64.b64encode(f.read()).decode("utf-8")
            scorecard_images[chain_name] = f"data:image/png;base64,{img_b64}"

    serialised = {k: {x: v for x, v in r.items() if x != "chain_obj"}
                 for k, r in results.items()}

    return jsonify({
        "success": True,
        "name": name,
        "results": serialised,
        "scorecard_images": scorecard_images,
    })


# ── Main ─────────────────────────────────────────────────────────────────────

if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(description="AbDev Flask API Server")
    parser.add_argument("--port", type=int, default=5001)
    parser.add_argument("--host", default="0.0.0.0")
    parser.add_argument("--debug", action="store_true")
    args = parser.parse_args()

    log.info(f"Starting AbDev API on {args.host}:{args.port}")
    app.run(host=args.host, port=args.port, debug=args.debug)
