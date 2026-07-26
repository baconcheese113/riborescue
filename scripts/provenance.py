"""Write the manifest a reader needs to regenerate this repository's results.

Reproducibility here rests on regeneration rather than stored bytes, so the manifest records what
would have to be re-fetched and re-run: the commit, the locked toolchain, every input with the
checksum it is verified against, the pinned upstream pipeline and container tags, the command that
produces each output, and each output's own hash.

Anything it cannot find is written into `missing` rather than omitted. A manifest that quietly
skips an absent input is a manifest that claims more than it checked.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import subprocess
import tomllib
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parent.parent

# Every fetched input, the release that identifies it, and the terms it comes under. Paths are
# relative to the repository root; each is hashed if present and reported as missing if not.
REFERENCES = [
    {
        "name": "GENCODE primary assembly genome",
        "release": "v50 (GRCh38)",
        "path": "data/gencode/GRCh38.primary_assembly.genome.fa.gz",
        "source": "https://www.gencodegenes.org/human/release_50.html",
        "license": "Open, per GENCODE/EMBL-EBI terms of use",
    },
    {
        "name": "GENCODE comprehensive annotation",
        "release": "v50",
        "path": "data/gencode/gencode.v50.primary_assembly.annotation.gtf.gz",
        "source": "https://www.gencodegenes.org/human/release_50.html",
        "license": "Open, per GENCODE/EMBL-EBI terms of use",
    },
    {
        "name": "GENCODE transcript sequences",
        "release": "v50",
        "path": "data/gencode/gencode.v50.transcripts.fa.gz",
        "source": "https://www.gencodegenes.org/human/release_50.html",
        "license": "Open, per GENCODE/EMBL-EBI terms of use",
    },
    {
        "name": "MANE Select genomic annotation",
        "release": "v1.5",
        "path": "data/mane/MANE.GRCh38.v1.5.refseq_genomic.gff.gz",
        "source": "https://ftp.ncbi.nlm.nih.gov/refseq/MANE/MANE_human/release_1.5/",
        "license": "US Government work, public domain",
    },
    {
        "name": "MANE Select transcripts",
        "release": "v1.5",
        "path": "data/mane/MANE.GRCh38.v1.5.refseq_rna.fna.gz",
        "source": "https://ftp.ncbi.nlm.nih.gov/refseq/MANE/MANE_human/release_1.5/",
        "license": "US Government work, public domain",
    },
    {
        "name": "MANE Select proteins",
        "release": "v1.5",
        "path": "data/mane/MANE.GRCh38.v1.5.refseq_protein.faa.gz",
        "source": "https://ftp.ncbi.nlm.nih.gov/refseq/MANE/MANE_human/release_1.5/",
        "license": "US Government work, public domain",
    },
    {
        "name": "MANE summary",
        "release": "v1.5",
        "path": "data/mane/MANE.GRCh38.v1.5.summary.txt.gz",
        "source": "https://ftp.ncbi.nlm.nih.gov/refseq/MANE/MANE_human/release_1.5/",
        "license": "US Government work, public domain",
    },
    {
        "name": "ClinVar variant summary",
        "release": "20260715",
        "path": "data/clinvar/clinvar_20260715.vcf.gz",
        "source": "https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh38/",
        "license": "US Government work, public domain",
    },
    {
        "name": "Human rDNA repeating unit",
        "release": "U13369.1",
        "path": "data/rdna/U13369.1.fasta",
        "source": "https://www.ncbi.nlm.nih.gov/nuccore/U13369.1",
        "license": "US Government work, public domain",
    },
    {
        "name": "Toledano et al. 2024 readthrough measurements",
        "release": "Extended Data Table 2, via the authors' Stop_codon_readthrough release",
        "path": "data/toledano/treated_samples.rds",
        "source": "https://github.com/lehner-lab/Stop_codon_readthrough",
        "license": "Per the publication and the authors' repository",
    },
]

# Each output the release stands behind, and the Pixi task that regenerates it.
OUTPUTS = [
    ("results/clinvar_nonsense.tsv", "clinvar"),
    ("results/clinvar_contexts.tsv", "contexts"),
    ("results/evaluation.tsv", "evaluate"),
    ("results/variant_therapy_scores.tsv", "score"),
    ("results/amenability_landscape.tsv", "landscape"),
    ("results/landscape_summary.tsv", "landscape"),
    ("results/nmd.tsv", "nmd"),
    ("results/aenmd_verdicts.tsv", "aenmd-verdicts"),
    ("results/nmdetective.tsv", "nmdetective"),
    ("results/diseases.tsv", "diseases"),
    ("results/disease_coverage.tsv", "disease-coverage"),
    ("results/trna_coverage.tsv", "trna-coverage"),
    ("results/trna_panel_genes.tsv", "trna-panel"),
    ("results/base_editing.tsv", "base-editing"),
    ("results/escape_map.tsv", "escape-map"),
    ("results/experiments.tsv", "experiments"),
    ("results/rnaseq/counts.tsv", "quantify"),
    ("results/rnaseq/tpm.tsv", "expression"),
    ("results/rnaseq/top_expressed.tsv", "expression"),
    ("results/rnaseq/top_expressed_nuclear.tsv", "expression"),
    ("results/rnaseq/composition.tsv", "expression"),
    ("results/staged_runs.tsv", "stage-runs"),
    ("results/reads/qc/multiqc_data/multiqc_general_stats.txt", "reads"),
    ("ledger/falsification.tsv", "hand-maintained"),
    ("experiments/programs.tsv", "hand-maintained"),
    ("frontend/public/riborescue.json", "export-web-safety"),
    ("frontend/public/riborescue_research.json", "export-research-tiers"),
]

# Containers the pipeline pins. Tags, not digests — resolving a digest needs a registry the release
# build does not require, so the gap is stated rather than papered over.
CONTAINERS = [
    "quay.io/biocontainers/fastqc:0.12.1--hdfd78af_0",
    "quay.io/biocontainers/cutadapt:5.2--py312hfabe715_2",
    "quay.io/biocontainers/star:2.7.11b--h5ca1c30_8",
    "quay.io/biocontainers/multiqc:1.35--pyhdfd78af_1",
    "ghcr.io/baconcheese113/riborescue:0.1.0",
]

# Dependencies a reader will not be able to install, and what the repository does instead.
UNAVAILABLE = [
    {
        "dependency": "predNMD",
        "why": "No public release; requested from the authors and not obtained",
        "consequence": "The NMD ensemble ships three predictors — the guideline rule, the full "
        "Lindeboom rule set, aenmd and NMDetective-AI — and no predNMD tier.",
    },
    {
        "dependency": "GSE179274 suppressor-tRNA contrasts",
        "why": "Five of ten libraries miss the pre-declared depth or P-site offset (ADR-0011, "
        "ADR-0013), including both suppressor-tRNA libraries and their EGFP control",
        "consequence": "The dataset is refused for confirmatory use. It is inconclusive under the "
        "frozen calibration gate, not negative. Only the exploratory detectability arm "
        "(ADR-0019) reads it.",
    },
    {
        "dependency": "Container image digests",
        "why": "The pipeline pins tags; digests are resolved by the registry at run time",
        "consequence": "Regeneration is reproducible to the tag, not to the byte. A digest pin is "
        "the remaining step.",
    },
]


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1 << 22), b""):
            digest.update(block)
    return digest.hexdigest()


def _git(*args: str) -> str:
    return subprocess.run(
        ["git", *args], cwd=ROOT, capture_output=True, text=True, check=True
    ).stdout.strip()


def _file_record(relative: str, missing: list[str], **extra: Any) -> dict[str, Any] | None:
    path = ROOT / relative
    if not path.exists():
        missing.append(relative)
        return None
    return {"path": relative, "bytes": path.stat().st_size, "sha256": _sha256(path), **extra}


def _sequencing_runs(missing: list[str]) -> list[dict[str, Any]]:
    """Every archived run, with the checksum the fetch verifies against.

    The published MD5 is the provenance statement, not a local re-hash: it is what makes a
    re-download the same bytes, and it holds whether or not the FASTQ is on this disk.
    """

    samples = ROOT / "pipeline/assets/riboseq_samples.tsv"
    if not samples.exists():
        missing.append(str(samples.relative_to(ROOT)))
        return []
    lines = samples.read_text().splitlines()
    header = lines[0].split("\t")
    runs = []
    for line in lines[1:]:
        row = dict(zip(header, line.split("\t"), strict=False))
        files = [
            {"url": row[url_key], "md5": row[md5_key]}
            for url_key, md5_key in (("fastq_1_url", "fastq_1_md5"), ("fastq_2_url", "fastq_2_md5"))
            if row.get(url_key)
        ]
        runs.append(
            {
                "dataset": row["dataset"],
                "sample": row["sample"],
                "run_accession": row["run_accession"],
                "assay": row["assay"],
                "layout": row["layout"],
                "cell_line": row["cell_line"],
                "treatment": row["treatment"],
                "read_count": int(row["read_count"]),
                "files": files,
            }
        )
    return runs


def _tasks() -> dict[str, str]:
    """The command behind each output, read from the manifest that actually runs it."""

    manifest = tomllib.loads((ROOT / "pyproject.toml").read_text())
    tasks = manifest["tool"]["pixi"]["feature"]["base"]["tasks"]
    out = {}
    for name, task in tasks.items():
        command = task if isinstance(task, str) else task.get("cmd", "")
        depends = [] if isinstance(task, str) else task.get("depends-on", [])
        out[name] = " ".join((command or f"depends-on: {', '.join(depends)}").split())
    return out


def build() -> dict[str, Any]:
    missing: list[str] = []
    tasks = _tasks()
    return {
        "repository": {
            "commit": _git("rev-parse", "HEAD"),
            "branch": _git("rev-parse", "--abbrev-ref", "HEAD"),
            "remote": _git("config", "--get", "remote.origin.url"),
            "uncommitted_changes": bool(_git("status", "--porcelain")),
        },
        "toolchain": {
            "platform": "linux-64",
            "manager": "pixi",
            "lockfile": _file_record("pixi.lock", missing),
            "manifest": _file_record("pyproject.toml", missing),
            "environments": ["default", "runtime", "psite", "aenmd", "nmdetective"],
        },
        "upstream_pipeline": {
            "pipeline": "nf-core/riboseq",
            "revision": "1.2.0",
            "consumed_as": "external pinned pipeline; this repository authors its own minimal "
            "workflow and no nf-core modules",
        },
        "containers": [{"image": image, "digest": None} for image in CONTAINERS],
        "references": [
            record
            for reference in REFERENCES
            if (
                record := _file_record(
                    reference["path"],
                    missing,
                    name=reference["name"],
                    release=reference["release"],
                    source=reference["source"],
                    license=reference["license"],
                )
            )
            is not None
        ],
        "sequencing_runs": _sequencing_runs(missing),
        "outputs": [
            record
            for path, task in OUTPUTS
            if (
                record := _file_record(
                    path, missing, task=task, command=tasks.get(task, "hand-maintained")
                )
            )
            is not None
        ],
        "unavailable": UNAVAILABLE,
        "missing": sorted(missing),
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out", type=Path, default=ROOT / "results/provenance.json")
    arguments = parser.parse_args()

    manifest = build()
    arguments.out.parent.mkdir(parents=True, exist_ok=True)
    arguments.out.write_text(json.dumps(manifest, indent=2) + "\n")

    print(f"{arguments.out.relative_to(ROOT)}")
    print(f"  commit          {manifest['repository']['commit'][:12]}")
    print(f"  references      {len(manifest['references'])} of {len(REFERENCES)}")
    print(f"  runs            {len(manifest['sequencing_runs'])}")
    print(f"  outputs         {len(manifest['outputs'])} of {len(OUTPUTS)}")
    if manifest["missing"]:
        print(f"  missing         {len(manifest['missing'])}")
        for path in manifest["missing"]:
            print(f"                  {path}")


if __name__ == "__main__":
    main()
