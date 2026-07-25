#!/usr/bin/env python
"""Score ClinVar nonsense variants with NMDetective-AI, on each gene's MANE Select transcript.

NMDetective-AI (Veiner et al. 2026; Vejni/NMDetectiveAI, MIT) is the model tier's deep member — a
fine-tuned Orthrus/Mamba encoder that predicts an NMD-efficiency score for a premature stop from the
transcript sequence. This encodes each variant on its MANE transcript into the model's six-track
input and runs it on the GPU, keeping every variant the model could not place (a transcript its
GENCODE v26 reference lacks, or an encoding error) explicitly unavailable with the reason.

Runs only in the `nmdetective` Pixi environment. Resumable: an existing output is read back and only
the variants not already scored are run, so an interrupted full pass continues where it stopped.
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import pandas as pd
import torch


def _mane_ensembl_by_gene(summary_path: str) -> dict[int, str]:
    summary = pd.read_csv(summary_path, sep="\t", dtype=str)
    gene = summary["#NCBI_GeneID"].str.removeprefix("GeneID:").astype(int)
    ensembl = summary["Ensembl_nuc"].str.split(".").str[0]
    select = summary["MANE_status"].str.contains("MANE Select", na=False)
    return dict(zip(gene[select], ensembl[select], strict=True))


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--nonsense", required=True)
    parser.add_argument("--mane-summary", required=True)
    parser.add_argument("--model", required=True)
    parser.add_argument("--out", required=True)
    parser.add_argument("--limit", type=int, default=None, help="Score only the first N (smoke).")
    parser.add_argument("--batch-log", type=int, default=500)
    args = parser.parse_args()

    import NMD.data.transcripts as transcripts
    from NMD.data.transcripts import create_six_track_encoding_with_variant
    from NMD.modeling.predict import _predict_batch, _setup_model
    from NMD.modeling.TrainerConfig import TrainerConfig
    from NMD.utils import load_model

    # The encoding reconstructs the GENCODE genome on every call; caching it by version leaves the
    # encoding math untouched but stops the reload from dominating the per-variant cost.
    _genome_cache: dict[str, object] = {}
    _load_genome = transcripts.Genome

    def _cached_genome(version: str) -> object:
        if version not in _genome_cache:
            _genome_cache[version] = _load_genome(version)
        return _genome_cache[version]

    transcripts.Genome = _cached_genome

    config = TrainerConfig()
    model, _, _, _, device = _setup_model(config)
    load_model(model, args.model, device=device)
    model.eval()
    print(f"model on {device}, cuda={torch.cuda.is_available()}", flush=True)

    mane = _mane_ensembl_by_gene(args.mane_summary)
    variants = pd.read_csv(args.nonsense, sep="\t", dtype={"chrom": str}).drop_duplicates(
        "variant_id"
    )
    if args.limit:
        variants = variants.head(args.limit)

    done: dict[str, dict] = {}
    out_path = Path(args.out)
    if out_path.exists():
        prior = pd.read_csv(out_path, sep="\t")
        done = {row.variant_id: row._asdict() for row in prior.itertuples(index=False)}
        print(f"resuming: {len(done)} variants already scored", flush=True)

    records: list[dict] = list(done.values())
    scored_this_run = 0
    for i, row in enumerate(variants.itertuples(index=False)):
        if row.variant_id in done:
            continue
        enst = mane.get(int(row.gene_id))
        base = {"variant_id": row.variant_id}
        if enst is None:
            records.append(base | {"nmd_available": False, "reason": "no_mane_ensembl"})
            continue
        var_str = f"{row.chrom}:{row.pos}:{row.ref}:{row.alt}"
        try:
            six = create_six_track_encoding_with_variant(enst, var_str, "gencode.v26")
            efficiency = float(_predict_batch(model, [six], device)[0])
            records.append(
                base | {"nmd_available": True, "reason": "", "nmd_efficiency": efficiency}
            )
        except Exception as error:
            # Any placement failure — a transcript GENCODE v26 lacks, an unplaceable variant — is an
            # honest unavailable, not a crash; the reason names the exception type.
            reason = type(error).__name__
            records.append(base | {"nmd_available": False, "reason": f"encode_error:{reason}"})
        scored_this_run += 1
        if scored_this_run % args.batch_log == 0:
            pd.DataFrame(records).to_csv(out_path, sep="\t", index=False)
            print(f"  scored {scored_this_run} this run ({i + 1}/{len(variants)})", flush=True)

    result = pd.DataFrame(records)
    result.to_csv(out_path, sep="\t", index=False)
    avail = int(result["nmd_available"].sum()) if "nmd_available" in result else 0
    print(f"wrote {len(result)} variants, {avail} scored → {out_path}", flush=True)
    sys.exit(0)


if __name__ == "__main__":
    main()
