"""Write the synthetic libraries the pipeline's read tests run over.

Two library shapes, because the trimming step treats them differently. A footprint library's insert
is shorter than the read, so every read runs into the linker; a transcriptome library's fragments
are mostly longer, so most reads never reach it. Both are needed to show that the adapter check
applies to the first and not the second.

    pixi run python scripts/make_read_fixtures.py
"""

import gzip
import random
from pathlib import Path

SEED = 721
READS = 200
RPF_ADAPTER = "CTGTAGGCACCATCAAT"
TRUSEQ_R1 = "AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC"
TRUSEQ_R2 = "AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTA"
FOOTPRINT_LENGTHS = (28, 29, 30, 31, 32)
FRAGMENT_LENGTHS = (60, 120, 120, 120)

DATA = Path(__file__).parents[1] / "pipeline/tests/data/reads"


def bases(rng: random.Random, count: int) -> str:
    return "".join(rng.choice("ACGT") for _ in range(count))


def write(path: Path, reads: list[tuple[str, str]]) -> None:
    with gzip.open(path, "wt", newline="\n") as handle:
        for name, read in reads:
            handle.write(f"@{name}\n{read}\n+\n{'I' * len(read)}\n")


def main() -> None:
    rng = random.Random(SEED)
    DATA.mkdir(parents=True, exist_ok=True)

    footprints = [
        (
            f"rpf_{i}",
            (bases(rng, rng.choice(FOOTPRINT_LENGTHS)) + RPF_ADAPTER + bases(rng, 50))[:50],
        )
        for i in range(READS)
    ]
    write(DATA / "riboseq_test_1.fastq.gz", footprints)

    forward, reverse = [], []
    for i in range(READS):
        insert = rng.choice(FRAGMENT_LENGTHS)
        forward.append((f"rna_{i}/1", (bases(rng, insert) + TRUSEQ_R1 + bases(rng, 100))[:100]))
        reverse.append((f"rna_{i}/2", (bases(rng, insert) + TRUSEQ_R2 + bases(rng, 100))[:100]))
    write(DATA / "rnaseq_test_1.fastq.gz", forward)
    write(DATA / "rnaseq_test_2.fastq.gz", reverse)


if __name__ == "__main__":
    main()
