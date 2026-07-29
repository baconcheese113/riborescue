/** The two arms of the pipeline and what each tool is there to do. One lane never touches raw reads;
 * the other processes ~18 GB of them. Naming both is what lets a reader judge the claims. */
export type Tool = { name: string; role: string };
export type Lane = { key: string; title: string; scale: string; intro: string; tools: Tool[] };

export const LANES: Lane[] = [
  {
    key: "measurement",
    title: "The measurement line",
    scale: "~18 GB of raw sequencing",
    intro:
      "Public ribosome-footprint and matched RNA sequencing, taken from raw reads through to positions on the genome resolved to a single codon.",
    tools: [
      { name: "SRA toolkit", role: "stages the public sequencing runs, verified against their checksums" },
      { name: "FastQC", role: "reads the raw quality: base scores, duplication, adapter contamination" },
      { name: "cutadapt", role: "trims sequencing adapters and drops reads outside the usable length range" },
      { name: "STAR", role: "aligns short footprints to the human genome across splice junctions" },
      { name: "samtools", role: "sorts, indexes and filters the alignments" },
      { name: "riboWaltz", role: "finds each footprint's reading position and checks the three-nucleotide rhythm" },
      { name: "featureCounts", role: "counts transcripts in the matched RNA data, for baseline abundance" },
      { name: "MultiQC", role: "collects every quality metric above into one report" },
      { name: "pigz", role: "decompresses the read archive in parallel" },
    ],
  },
  {
    key: "prediction",
    title: "The prediction line",
    scale: "70,376 scored variants",
    intro:
      "Variant databases and genomic sequence context. This lane never touches a raw read, and runs on a laptop.",
    tools: [
      { name: "cyvcf2", role: "parses the ClinVar variant file" },
      { name: "Pydantic + Pandera", role: "enforces the data contracts that stop coordinates drifting between modules" },
      { name: "Biopython", role: "codon tables, translation and the sequence vocabulary" },
      { name: "bioframe", role: "overlaps genomic intervals — chosen over the alternative after it needed 13 GiB" },
      { name: "statsmodels + SciPy", role: "fits the readthrough model and runs the resampling behind every interval" },
      { name: "aenmd", role: "a published rule engine for whether the transcript survives" },
      { name: "NMDetective-AI", role: "a deep model on the same question, run on a GPU as an independent third opinion" },
      { name: "BLOSUM62", role: "scores how chemically similar the inserted amino acid is to the original" },
      { name: "Editor panel", role: "BE4max and ABE7.10 geometry against SpCas9 targeting rules" },
    ],
  },
  {
    key: "shared",
    title: "Underneath both",
    scale: "one lockfile",
    intro:
      "The parts that make a number reproducible rather than merely produced.",
    tools: [
      { name: "Nextflow", role: "orders the steps and resumes them; the scientific decisions stay in Python" },
      { name: "Pixi", role: "one lockfile across Python, R, the JVM and the command-line tools" },
      { name: "R + caret", role: "the independent implementation the Python model is checked against" },
      { name: "matplotlib", role: "the report figures, drawn from the same tables and palette as this site" },
    ],
  },
];
