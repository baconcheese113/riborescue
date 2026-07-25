"""The riborescue command line — every scientific step the pipeline runs enters through here."""

import subprocess
from pathlib import Path

import click
import pandas as pd
import pandera.errors
import pandera.pandas

from riborescue._version import __version__
from riborescue.core.contracts import Consequence, EvalConfig
from riborescue.core.handoff import UpstreamHandoff
from riborescue.core.inputs import INPUTS, UnknownInputError, data_root, fetch
from riborescue.core.tables import (
    PathogenicNonsense,
    ReadthroughLabels,
    SequencingRuns,
    StagedRuns,
    TriageInput,
    TriageOutput,
    read_table,
    write_table,
)
from riborescue.riboseq.calibration import read_manifest, select_lengths
from riborescue.riboseq.contaminants import write_contaminants
from riborescue.riboseq.expression import (
    composition,
    gene_symbols,
    library_depth,
    read_counts,
    top_expressed,
    tpm,
)
from riborescue.riboseq.native_stop_atlas import native_stop_occupancy, translate_extension
from riborescue.riboseq.reads import (
    ADAPTER_REACHED_BY,
    AdapterNotFoundError,
    summarise_alignment,
    summarise_trimming,
    survey_adapters,
)
from riborescue.riboseq.readthrough_assay import (
    PROGRAMMED_READTHROUGH,
    collapse_lengths,
    extension_windows,
    library_ratios,
    mane_select_transcripts,
    overlapping_downstream_cds,
    paired_effect,
    qualifying,
    signature,
    stalling,
    transcript_genes,
    unpaired_effect,
)
from riborescue.riboseq.sequencing import FASTQ_SUBDIR, stage
from riborescue.variants.aenmd import (
    aenmd_verdicts,
    mane_ensembl_by_gene,
    model_agreement,
    read_aenmd_rules,
)
from riborescue.variants.clinvar import pathogenic_nonsense
from riborescue.variants.context import contexts_for, disagreements_with_protein
from riborescue.variants.disease_coverage import disease_coverage, disease_reach_frontier
from riborescue.variants.diseases import normalize_conditions
from riborescue.variants.evaluation import (
    SEED,
    UnsupportedEvalConfigError,
    bootstrap_ci,
    evaluate,
    split,
)
from riborescue.variants.landscape import TOLERABLE_SHARE, Thresholds, landscape, summarise
from riborescue.variants.native_stop_predictions import (
    concordance,
    four_quadrants,
    native_stop_features,
)
from riborescue.variants.nmd_rules import disagreement_atlas, nmd_predictors
from riborescue.variants.readthrough_model import fit, relative_error_quantile
from riborescue.variants.research_export import build_research_aggregate
from riborescue.variants.residue import coverage_by_design
from riborescue.variants.suppressor_panels import coverage_frontier
from riborescue.variants.transcripts import load_transcripts, read_sequences
from riborescue.variants.triage import classify, classify_table
from riborescue.variants.web_export import (
    build_web_table,
    diverse_sample,
    read_web_inputs,
    safety_summary,
)

__all__ = ["main"]

_IN = click.argument("table", type=click.Path(exists=True, dir_okay=False, path_type=Path))
_OUT = click.option(
    "--out",
    required=True,
    type=click.Path(dir_okay=False, writable=True, path_type=Path),
    help="Where to write the resulting table.",
)


@click.group()
@click.version_option(__version__)
def main() -> None:
    """Match nonsense variants to candidate readthrough therapies."""


@main.command()
@click.argument("consequence", type=click.Choice([c.value for c in Consequence]))
@click.option(
    "--transcript-supported/--transcript-unsupported",
    default=True,
    help="Whether the variant's transcript is in the supported set.",
)
def triage(consequence: str, transcript_supported: bool) -> None:
    """Say whether readthrough applies to a variant CONSEQUENCE."""

    result = classify(Consequence(consequence), transcript_supported=transcript_supported)
    verdict = "applies" if result.applies else "does not apply"
    click.echo(f"{result.triage_class.value}: readthrough {verdict} — {result.reason}")


@main.command("triage-table")
@_IN
@_OUT
def triage_table(table: Path, out: Path) -> None:
    """Triage every variant in TABLE, writing the verdicts alongside them."""

    variants = _validated(TriageInput, read_table(table), table)
    triaged = _validated(TriageOutput, classify_table(variants), table)
    write_table(triaged, out)
    applies = int(triaged["applies"].sum())
    click.echo(f"triaged {len(triaged)} variants; readthrough applies to {applies}")


@main.command("validate-labels")
@_IN
def validate_labels(table: Path) -> None:
    """Check that TABLE holds well-formed readthrough efficiency labels."""

    labels = _validated(ReadthroughLabels, read_table(table), table)
    censored = int(labels["censored"].sum())
    click.echo(f"{table}: {len(labels)} labels, {censored} at the assay ceiling")


@main.command("validate-handoff")
@click.argument("manifest", type=click.Path(exists=True, dir_okay=False, path_type=Path))
@click.option(
    "--check-files/--no-check-files",
    default=False,
    help="Also require every declared upstream output to be present on disk.",
)
@click.option(
    "--results-root",
    type=click.Path(file_okay=False, path_type=Path),
    help="Where the results tree sits now, overriding the root the manifest names.",
)
def validate_handoff(manifest: Path, check_files: bool, results_root: Path | None) -> None:
    """Check the upstream handoff MANIFEST naming the nf-core/riboseq outputs to consume."""

    try:
        handoff = UpstreamHandoff.from_json(manifest)
    except ValueError as error:
        raise click.ClickException(
            f"{manifest} is not a valid handoff manifest\n{error}"
        ) from error
    if results_root is not None:
        handoff = handoff.model_copy(update={"results_root": results_root})
    if check_files and (missing := handoff.missing()):
        listed = "\n".join(f"  {name}: {path}" for name, path in missing)
        raise click.ClickException(f"{manifest} declares outputs that are not present\n{listed}")
    declared = len(tuple(handoff.outputs()))
    click.echo(f"{manifest}: {handoff.pipeline} {handoff.revision}, {declared} declared outputs")


@main.command("fetch")
@click.argument("names", nargs=-1, type=click.Choice(list(INPUTS)))
@click.option(
    "--data-root",
    "data_root_",
    type=click.Path(file_okay=False, path_type=Path),
    help="Where inputs are kept, overriding RIBORESCUE_DATA.",
)
@click.option("--force", is_flag=True, help="Fetch again even if the file already verifies.")
def fetch_inputs(names: tuple[str, ...], data_root_: Path | None, force: bool) -> None:
    """Fetch declared public inputs by NAME, or all of them, and verify their checksums."""

    root = data_root_ if data_root_ is not None else data_root()
    for name in names or tuple(INPUTS):
        declared = INPUTS[name]
        click.echo(f"{name}: {declared.source} ({declared.licence})")
        try:
            path = fetch(name, root, force=force)
        except (UnknownInputError, OSError) as error:
            raise click.ClickException(str(error)) from error
        click.echo(f"  {path} verified")


@main.command("stage-runs")
@click.argument("samplesheet", type=click.Path(exists=True, dir_okay=False, path_type=Path))
@_OUT
@click.option(
    "--data-root",
    "data_root_",
    type=click.Path(file_okay=False, path_type=Path),
    help="Where inputs are kept, overriding RIBORESCUE_DATA.",
)
@click.option("--force", is_flag=True, help="Fetch again even if the file already verifies.")
def stage_runs(samplesheet: Path, out: Path, data_root_: Path | None, force: bool) -> None:
    """Fetch the FASTQ named by SAMPLESHEET and write a sheet naming them on disk."""

    runs = _validated(SequencingRuns, read_table(samplesheet), samplesheet)
    root = data_root_ if data_root_ is not None else data_root()
    try:
        staged = stage(runs, root, force=force)
    except OSError as error:
        raise click.ClickException(str(error)) from error
    write_table(_validated(StagedRuns, staged, out), out)
    click.echo(f"{len(staged)} runs staged under {root / FASTQ_SUBDIR}")


@main.command("select-lengths")
@click.argument("dataset")
@click.option(
    "--frames",
    required=True,
    type=click.Path(exists=True, dir_okay=False, path_type=Path),
    help="frame_by_length.tsv from the survey calibration pass.",
)
@click.option(
    "--offsets",
    required=True,
    type=click.Path(exists=True, dir_okay=False, path_type=Path),
    help="psite_offsets.tsv from the same pass.",
)
@click.option(
    "--script-md5",
    type=click.Path(exists=True, dir_okay=False, path_type=Path),
    help="The checksum file the calibration driver wrote beside its results.",
)
@_OUT
def select_calibration_lengths(
    dataset: str, frames: Path, offsets: Path, script_md5: Path | None, out: Path
) -> None:
    """Choose DATASET's footprint lengths and judge every library on them.

    The set is one shared choice for the whole dataset, made from periodicity alone and before any
    contrast is run. A library that fails a predeclared threshold is not dropped and the threshold
    is not moved: the manifest records the failure and the assay refuses to run.
    """

    checksum = script_md5.read_text().strip() if script_md5 is not None else None
    try:
        manifest = select_lengths(dataset, read_table(frames), read_table(offsets), checksum)
    except ValueError as error:
        raise click.ClickException(str(error)) from error

    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(manifest.to_json())

    kept = ", ".join(str(length) for length in manifest.lengths) or "none"
    click.echo(f"{dataset}: lengths {kept}")
    for library in manifest.libraries:
        state = "pass" if library.passes else "FAIL"
        click.echo(
            f"  {state} {library.sample}: {library.psites:,} P-sites, "
            f"frame-0 {library.frame0_share:.1%}, {library.dominant_length} nt dominant, "
            f"offset {library.offset_from_5} nt"
        )
        for failure in library.failures:
            click.echo(f"       {failure}")
    if not manifest.passes:
        raise click.ClickException(f"{dataset} is inconclusive; its thresholds are not to be moved")


@main.command("atlas")
@click.argument("counts", type=click.Path(exists=True, dir_okay=False, path_type=Path))
@click.option("--gtf", required=True, type=click.Path(exists=True, path_type=Path))
@click.option(
    "--samplesheet",
    required=True,
    type=click.Path(exists=True, dir_okay=False, path_type=Path),
    help="The staged runs, naming each library's dataset and treatment.",
)
@click.option("--dataset", required=True, help="The experiment whose libraries the atlas measures.")
@click.option(
    "--control", required=True, help="The untreated arm a native stop should be tight in."
)
@click.option(
    "--transcripts",
    type=click.Path(exists=True, dir_okay=False, path_type=Path),
    help="Transcript FASTA, to translate the peptide each readthrough would add.",
)
@click.option(
    "--annotation",
    type=click.Path(exists=True, dir_okay=False, path_type=Path),
    help="The riboWaltz annotation table, for the native stop's position.",
)
@click.option(
    "--extensions",
    type=click.Path(exists=True, dir_okay=False, path_type=Path),
    help="The extension windows, for how far each readthrough runs.",
)
@click.option("--lengths", default="28:35", show_default=True, help="Footprint lengths to pool.")
@_OUT
def atlas_command(
    counts: Path,
    gtf: Path,
    samplesheet: Path,
    dataset: str,
    control: str,
    transcripts: Path | None,
    annotation: Path | None,
    extensions: Path | None,
    lengths: str,
    out: Path,
) -> None:
    """Measure readthrough past every transcript's NATIVE stop, per treatment arm.

    The empirical layer of the safety atlas: the same window that measured a drug's effect past
    premature stops, read at the native stop of each transcript. It is an anchor, not ground truth —
    it covers only what the cell line expressed at depth, and its per-transcript counts are small.
    """

    lo, hi = (int(part) for part in lengths.split(":"))
    footprints = range(lo, hi + 1)

    runs = _validated(StagedRuns, read_table(samplesheet), samplesheet)
    library = runs.loc[
        (runs["assay"] == "riboseq") & (runs["dataset"] == dataset),
        ["sample", "treatment"],
    ]
    arms = {str(t): list(rows["sample"]) for t, rows in library.groupby("treatment")}
    if control not in arms:
        raise click.ClickException(f"{samplesheet} names no {control} arm for {dataset}")

    genes = transcript_genes(gtf)
    programmed = frozenset(genes.index[genes.isin(PROGRAMMED_READTHROUGH)])
    excluded = overlapping_downstream_cds(gtf) | programmed

    measured = read_table(counts)
    table = native_stop_occupancy(measured, arms, footprints, excluded_transcripts=excluded)
    table["gene"] = table["transcript"].map(genes)
    table["control_occupancy"] = table[f"{control}_occupancy"]

    if transcripts is not None and annotation is not None and extensions is not None:
        peptides = translate_extension(transcripts, read_table(annotation), read_table(extensions))
        table = table.merge(peptides, on="transcript", how="left")

    write_table(table, out)
    included = table[table["included"]]
    tight = included[included["control_occupancy"] < 0.05]
    click.echo(
        f"{len(table):,} transcripts, {len(included):,} deep enough, "
        f"{len(tight):,} with a native stop tight in {control}"
    )
    for arm in sorted(arms):
        if arm == control:
            continue
        lift = (included[f"{arm}_occupancy"] - included["control_occupancy"]).mean()
        click.echo(f"  {arm}: mean downstream occupancy {lift:+.4f} over {control}")


@main.command("atlas-predict")
@click.option(
    "--transcripts",
    required=True,
    type=click.Path(exists=True, dir_okay=False, path_type=Path),
    help="Transcript FASTA.",
)
@click.option(
    "--annotation",
    required=True,
    type=click.Path(exists=True, dir_okay=False, path_type=Path),
    help="The riboWaltz annotation table, for each native stop's position.",
)
@click.option(
    "--training",
    required=True,
    type=click.Path(exists=True, dir_okay=False, path_type=Path),
    help="The G418 training features the amenability model is fit on.",
)
@click.option(
    "--measured",
    required=True,
    type=click.Path(exists=True, dir_okay=False, path_type=Path),
    help="The empirical atlas table, for the measured occupancy to compare against.",
)
@click.option(
    "--gtf",
    required=True,
    type=click.Path(exists=True, path_type=Path),
    help="The GENCODE annotation, to pick the MANE Select transcript per gene.",
)
@_OUT
def atlas_predict(
    transcripts: Path, annotation: Path, training: Path, measured: Path, gtf: Path, out: Path
) -> None:
    """Score every canonical stop with the model, and rank-compare it to the measured occupancy.

    The predicted layer of the safety atlas, and an extrapolation by design: the model was fit on
    premature stops, so applying it to a canonical stop is out of distribution even where the stop's
    own feature levels were seen in training. The prediction is kept in its own column and never
    blended with the measurement. A prediction row is one GENCODE coding transcript; the concordance
    is computed only over MANE Select transcripts, one canonical stop per gene, so a gene's isoforms
    do not pseudo-replicate a shared context.
    """

    features = native_stop_features(transcripts, read_table(annotation)).set_index("transcript")
    model = fit(read_table(training))
    supported = model.known_levels(features)
    predicted = model.predict(features[supported])

    table = pd.DataFrame(
        {
            "predicted_g418": predicted.reindex(features.index),
            "application_domain": "native_stop_extrapolation",
            "feature_support": ["supported" if k else "unsupported" for k in supported],
        },
        index=features.index,
    ).reset_index()

    atlas = read_table(measured)
    atlas = atlas[atlas["included"]].copy()
    atlas["measured_lift"] = atlas["g418_occupancy"] - atlas["control_occupancy"]
    combined = table.merge(
        atlas[
            [
                "transcript",
                "gene",
                "control_occupancy",
                "g418_occupancy",
                "g418_depth",
                "measured_lift",
            ]
        ],
        on="transcript",
        how="left",
    )
    combined["measurable"] = combined["measured_lift"].notna()
    combined["mane_select"] = combined["transcript"].isin(mane_select_transcripts(gtf))

    # One canonical stop per gene: isoforms sharing a stop context must not each count as evidence.
    primary = combined[combined["mane_select"]].dropna(subset=["predicted_g418", "measured_lift"])
    stats = concordance(primary["predicted_g418"], primary["measured_lift"])

    # Thresholds are chosen here, after the fact, so the four groups describe the data rather than
    # classify against a pre-registered line: the predicted cut is the upper quartile, the measured
    # cut a five-point occupancy lift.
    cut_p = float(primary["predicted_g418"].quantile(0.75)) if len(primary) else 0.0
    cut_m = 0.05
    combined["group"] = ""
    combined.loc[primary.index, "group"] = four_quadrants(
        primary["predicted_g418"], primary["measured_lift"], cut_p, cut_m
    )

    out.parent.mkdir(parents=True, exist_ok=True)
    write_table(combined, out)

    click.echo(
        f"{len(table):,} canonical stops scored, one per GENCODE coding transcript "
        f"({int(supported.sum()):,} feature-supported, {int((~supported).sum()):,} unsupported)"
    )
    click.echo(
        f"rank concordance over {stats['n']:,} MANE Select stops (one per gene): "
        f"Spearman {stats['rho']} [{stats['low']}, {stats['high']}]"
    )
    click.echo(f"  four groups at predicted≥{cut_p:.3f}, measured lift≥{cut_m:.2f} (descriptive):")
    counts = combined.loc[primary.index, "group"].value_counts()
    for group in ("both", "predicted only", "measured only", "neither"):
        click.echo(f"    {group:>14}: {int(counts.get(group, 0)):,}")


@main.command("export-web")
@click.option(
    "--landscape",
    required=True,
    type=click.Path(exists=True, dir_okay=False, path_type=Path),
    help="The per-variant landscape table.",
)
@click.option(
    "--amenability",
    required=True,
    type=click.Path(exists=True, dir_okay=False, path_type=Path),
    help="The per-therapy amenability scores.",
)
@click.option(
    "--predicted",
    type=click.Path(exists=True, dir_okay=False, path_type=Path),
    help="The native-stop prediction table, to summarise the safety atlas for the viewer.",
)
@click.option(
    "--nmd",
    "nmd_table",
    type=click.Path(exists=True, dir_okay=False, path_type=Path),
    help="The NMD predictors table, to show each variant's two escape verdicts and their rules.",
)
@click.option(
    "--aenmd",
    "aenmd_table",
    type=click.Path(exists=True, dir_okay=False, path_type=Path),
    help="The aenmd verdicts table, to show the published tool's call beside the rule tier.",
)
@click.option(
    "--nmdetective",
    "nmdetective_table",
    type=click.Path(exists=True, dir_okay=False, path_type=Path),
    help="The NMDetective-AI table, to show the deep model's efficiency score per variant.",
)
@click.option(
    "--sample",
    type=int,
    help="Write a diverse sample of this many variants — the example artifact — instead of all.",
)
@_OUT
def export_web(
    landscape: Path,
    amenability: Path,
    predicted: Path | None,
    nmd_table: Path | None,
    aenmd_table: Path | None,
    nmdetective_table: Path | None,
    sample: int | None,
    out: Path,
) -> None:
    """Build the compact JSON the web app reads from the pipeline's result tables."""

    land, amen = read_web_inputs(landscape, amenability)
    variant_ids = diverse_sample(land, amen, sample) if sample is not None else None
    therapies = sorted(amen["therapy_id"].unique())
    safety = safety_summary(read_table(predicted), therapies) if predicted is not None else None
    nmd = read_table(nmd_table) if nmd_table is not None else None
    aenmd = read_table(aenmd_table) if aenmd_table is not None else None
    nmdetective = read_table(nmdetective_table) if nmdetective_table is not None else None
    table = build_web_table(
        land,
        amen,
        variant_ids=variant_ids,
        safety=safety,
        nmd=nmd,
        aenmd=aenmd,
        nmdetective=nmdetective,
    )
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(table.to_json())
    click.echo(
        f"{len(table.variants):,} variants, {len(table.therapies)} therapies"
        + (f", safety atlas over {safety['analysed']:,} stops" if safety else "")
        + f" → {out}"
    )


@main.command("adapter-survey")
@click.argument("samplesheet", type=click.Path(exists=True, dir_okay=False, path_type=Path))
@_OUT
def adapter_survey(samplesheet: Path, out: Path) -> None:
    """Check every footprint library in SAMPLESHEET against the adapter it declares.

    Run on the staged reads, before trimming. What a series says about its adapter has been wrong
    or absent three times out of four, and trimming the wrong sequence is silent: the reads align,
    softly clipped, carrying linker into every P-site.
    """

    runs = _validated(StagedRuns, read_table(samplesheet), samplesheet)
    try:
        survey = survey_adapters(runs)
    except (AdapterNotFoundError, OSError, ValueError) as error:
        raise click.ClickException(str(error)) from error
    write_table(survey, out)
    for row in survey.itertuples():
        click.echo(
            f"{row.sample}: adapter in {row.adapter_rate:.1%} of reads, "
            f"footprint {row.footprint_p10}-{row.footprint_p90} nt "
            f"(median {row.footprint_median})"
        )


@main.command("trim-summary")
@click.argument(
    "reports", nargs=-1, required=True, type=click.Path(exists=True, dir_okay=False, path_type=Path)
)
@click.option(
    "--samplesheet",
    required=True,
    type=click.Path(exists=True, dir_okay=False, path_type=Path),
    help="The staged runs the reports came from, naming each library's assay.",
)
@_OUT
def trim_summary(reports: tuple[Path, ...], samplesheet: Path, out: Path) -> None:
    """Summarise cutadapt REPORTS into raw-against-cleaned counts per library."""

    runs = _validated(StagedRuns, read_table(samplesheet), samplesheet)
    expected = set(runs.loc[runs["assay"].isin(ADAPTER_REACHED_BY), "sample"])
    try:
        summary = summarise_trimming(sorted(reports), expected)
    except AdapterNotFoundError as error:
        raise click.ClickException(str(error)) from error
    write_table(summary, out)
    for row in summary.itertuples():
        click.echo(
            f"{row.sample}: {row.reads_raw:,} raw, {row.reads_cleaned:,} cleaned "
            f"({row.reads_retained:.1%} retained, adapter in {row.adapter_rate:.1%})"
        )


@main.command("contaminants")
@click.argument("transcripts", type=click.Path(exists=True, dir_okay=False, path_type=Path))
@click.option(
    "--include",
    multiple=True,
    type=click.Path(exists=True, dir_okay=False, path_type=Path),
    help="Further FASTA to deplete against, such as the rDNA repeating unit.",
)
@_OUT
def contaminant_sequences(transcripts: Path, include: tuple[Path, ...], out: Path) -> None:
    """Write the structural RNA in a GENCODE TRANSCRIPTS FASTA, the sequences to deplete."""

    try:
        written = write_contaminants(transcripts, out, include)
    except (OSError, ValueError) as error:
        raise click.ClickException(str(error)) from error
    click.echo(f"{written:,} contaminant sequences written to {out}")


@main.command("alignment-summary")
@click.argument(
    "logs", nargs=-1, required=True, type=click.Path(exists=True, dir_okay=False, path_type=Path)
)
@_OUT
def alignment_summary(logs: tuple[Path, ...], out: Path) -> None:
    """Summarise STAR's final LOGS into per-library alignment metrics."""

    summary = summarise_alignment(sorted(logs))
    write_table(summary, out)
    for row in summary.itertuples():
        click.echo(
            f"{row.sample}: {row.reads_input:,.0f} reads, {row.unique_rate:.1f}% unique, "
            f"{row.multimapped_rate:.1f}% multimapped, "
            f"mean mapped length {row.mapped_length_mean:.1f}"
        )


@main.command("expression")
@click.argument("counts", type=click.Path(exists=True, dir_okay=False, path_type=Path))
@click.option(
    "--annotation",
    required=True,
    type=click.Path(exists=True, dir_okay=False, path_type=Path),
    help="The GTF the counts were made against, for gene symbols.",
)
@click.option(
    "--top", default=20, show_default=True, help="How many of the most expressed genes to write."
)
@click.option(
    "--tpm-out",
    required=True,
    type=click.Path(dir_okay=False, writable=True, path_type=Path),
    help="Where to write the per-gene TPM table.",
)
@click.option(
    "--nuclear-out",
    type=click.Path(dir_okay=False, writable=True, path_type=Path),
    help="Where to write the same ranking with mitochondrial and structural RNA set aside.",
)
@click.option(
    "--composition-out",
    type=click.Path(dir_okay=False, writable=True, path_type=Path),
    help="Where to write each library's mitochondrial, structural and nuclear shares.",
)
@_OUT
def expression(
    counts: Path,
    annotation: Path,
    top: int,
    tpm_out: Path,
    nuclear_out: Path | None,
    composition_out: Path | None,
    out: Path,
) -> None:
    """Turn featureCounts COUNTS from the matched RNA-seq arm into TPM and the most expressed genes.

    Abundance, not occupancy: this is the RNA-seq arm, and it answers how much message there is
    rather than where ribosomes sit. TPM is normalised per library, so a library carrying more
    surviving mitochondrial or structural RNA deflates every other gene in it — which is why the
    composition is reported beside the ranking, and why the ranking is offered both ways.
    """

    table = read_counts(counts)
    depth = library_depth(table)
    expression_table = tpm(table)
    symbols = gene_symbols(annotation)
    ranked = top_expressed(expression_table, symbols, top)

    written = [(tpm_out, expression_table), (out, ranked)]
    if nuclear_out is not None:
        written.append(
            (
                nuclear_out,
                top_expressed(expression_table, symbols, top, ("mitochondrial", "structural")),
            )
        )
    shares = composition(expression_table, symbols)
    if composition_out is not None:
        written.append((composition_out, shares))
    for path, frame in written:
        path.parent.mkdir(parents=True, exist_ok=True)
        write_table(frame, path)

    click.echo(f"{len(table):,} genes over {len(depth)} libraries")
    for library, assigned in depth.items():
        click.echo(f"  {library}: {assigned:,} assigned")
    click.echo("library composition (% of TPM):")
    libraries = [column for column in shares.columns if column != "gene_class"]
    for row in shares.itertuples(index=False):
        values = "  ".join(f"{getattr(row, library):>5.1f}" for library in libraries)
        click.echo(f"  {row.gene_class:<14} {values}")
    click.echo(f"top {len(ranked)} by mean TPM:")
    for row in ranked.head(5).itertuples():
        click.echo(f"  {row.gene_symbol:<12} {row.mean_tpm:>12,.0f}")


@main.command("evaluate")
@click.argument(
    "features",
    nargs=-1,
    required=True,
    type=click.Path(exists=True, dir_okay=False, path_type=Path),
)
@click.option(
    "--config",
    "configs",
    multiple=True,
    type=click.Choice([c.value for c in EvalConfig]),
    default=[EvalConfig.published_random_cv.value, EvalConfig.grouped_by_gene.value],
    help="Which named evaluation protocols to run.",
)
@click.option("--seed", default=SEED, show_default=True, help="Seed the splits are drawn from.")
@click.option("--out", type=click.Path(dir_okay=False, path_type=Path), help="Write rounds here.")
def evaluate_features(
    features: tuple[Path, ...], configs: tuple[str, ...], seed: int, out: Path | None
) -> None:
    """Score the baseline on each FEATURES table under each named evaluation protocol.

    A dataset is named after its file, so `features_G418.tsv.gz` reports as G418.
    """

    rounds = []
    # Sorted so the output does not depend on the order the caller's shell or workflow globbed them.
    for path in sorted(features):
        table = read_table(path).set_index("row")
        dataset = path.name.split(".")[0].removeprefix("features_")
        for name in configs:
            config = EvalConfig(name)
            try:
                scored = evaluate(table, split(table, config, seed=seed))
            except UnsupportedEvalConfigError as error:
                raise click.ClickException(str(error)) from error
            interval = bootstrap_ci(scored["r2"], seed=seed)
            rounds.append(scored.assign(dataset=dataset, config=config.value))
            click.echo(
                f"{dataset:<12} {config.value:<22} "
                f"r² {interval.point:.3f}  95% CI [{interval.low:.3f}, {interval.high:.3f}]"
            )
    if out is not None:
        out.parent.mkdir(parents=True, exist_ok=True)
        write_table(pd.concat(rounds), out)
        click.echo(f"wrote {out}")


@main.command("clinvar")
@click.argument("vcf", type=click.Path(exists=True, dir_okay=False, path_type=Path))
@_OUT
def clinvar_variants(vcf: Path, out: Path) -> None:
    """Extract the pathogenic nonsense substitutions from a ClinVar VCF."""

    found = pathogenic_nonsense(vcf)
    variants = _validated(PathogenicNonsense, found.variants, vcf)
    out.parent.mkdir(parents=True, exist_ok=True)
    write_table(variants, out)
    genes = variants["gene_symbol"].nunique()
    reviewed = int((variants["review_stars"] >= 2).sum())
    click.echo(
        f"{len(variants)} pathogenic nonsense variants across {genes} genes; "
        f"{reviewed} with two or more review stars"
    )
    if found.ambiguous_alleles:
        click.echo(f"excluded {found.ambiguous_alleles} with an ambiguous alternate allele")


@main.command("diseases")
@_IN
@_OUT
def diseases(table: Path, out: Path) -> None:
    """Normalize each variant's ClinVar conditions to MedGen/OMIM/Orphanet identifiers.

    One row per variant-condition, keyed on MedGen. Placeholder and partially-mapped conditions are
    kept and labelled rather than dropped, so a denominator excludes them deliberately.
    """

    variants = read_table(table)
    normalized = normalize_conditions(variants)
    out.parent.mkdir(parents=True, exist_ok=True)
    write_table(normalized, out)

    real = normalized[~normalized["mapping_status"].isin(["placeholder", "unmapped"])]
    click.echo(
        f"{len(normalized)} variant-condition rows over "
        f"{variants['variant_id'].nunique()} variants; "
        f"{real['medgen'].nunique()} distinct MedGen condition entities (not all diseases)"
    )
    for status, count in normalized["mapping_status"].value_counts().items():
        click.echo(f"  {status:<12} {count}")


_CONTEXTS = click.option(
    "--contexts",
    required=True,
    type=click.Path(exists=True, dir_okay=False, path_type=Path),
    help="The contexts table, for which variants are scoreable and at which stop and residue.",
)


@main.command("disease-coverage")
@_IN
@_CONTEXTS
@_OUT
def disease_coverage_cmd(table: Path, contexts: Path, out: Path) -> None:
    """Per-condition-entity model coverage over each entity's eligible nonsense-variant denominator.

    A MedGen concept is a ClinVar condition — a disease, but possibly a finding or susceptibility —
    so these are condition entities, not verified diseases. Three metrics are kept apart: reach, the
    covered fraction, and complete coverage. Makes no claim about unmet therapeutic need.
    """

    coverage = disease_coverage(read_table(table), read_table(contexts))
    out.parent.mkdir(parents=True, exist_ok=True)
    write_table(coverage, out)

    reached = int(coverage["reach"].sum())
    complete = int(coverage["complete"].sum())
    click.echo(
        f"{len(coverage)} condition entities (MedGen concepts, not all diseases); "
        f"{reached} reached by ≥1 design, {complete} with every eligible variant covered"
    )
    head = coverage.head(5)
    for i in range(len(head)):
        row = head.iloc[i]
        click.echo(
            f"  {str(row['condition_name'])[:40]:<40} "
            f"{int(row['model_covered']):>4}/{int(row['eligible_variants']):<4} "
            f"({row['covered_fraction'] * 100:>5.1f}%)"
        )


@main.command("disease-panel")
@_IN
@_CONTEXTS
@_OUT
def disease_panel_cmd(table: Path, contexts: Path, out: Path) -> None:
    """Greedy panel reaching the most condition entities — one with any model-covered variant.

    This is the reach frontier, distinct from the per-entity coverage fraction: a design reaches an
    entity when it restores at least one of its variants. Reaching every entity is partly a closure
    property of the design universe, not a therapeutic result.
    """

    diseases = read_table(table)
    frontier = disease_reach_frontier(diseases, read_table(contexts))
    out.parent.mkdir(parents=True, exist_ok=True)
    write_table(frontier, out)

    eligible = diseases[diseases["mapping_status"].isin(["mapped", "medgen_only"])]
    entities = int(eligible["medgen"].nunique())
    pairs = len(eligible)
    excluded = int(len(diseases) - pairs)
    total = int(frontier["cumulative"].iloc[-1]) if len(frontier) else 0

    click.echo(
        f"{len(frontier)} designs reach every one of {total} condition entities represented by an "
        f"exact-restorable ClinVar nonsense SNV"
    )
    click.echo(
        f"  denominator: {entities} eligible condition entities over {pairs} unique "
        f"variant-condition pairs; {excluded} placeholder/unmapped rows excluded, "
        f"one-to-many deduplicated to distinct MedGen; {entities - total} entities have no "
        f"restorable variant"
    )
    click.echo("  reaching all is partly a closure property of the design universe, not a result")
    for k in (1, 3, 5, 10):
        if k <= len(frontier):
            row = frontier.iloc[k - 1]
            click.echo(
                f"  panel of {k:>2}: {int(row['cumulative']):>4} of {total} entities "
                f"({row['cumulative_fraction'] * 100:>5.1f}%), "
                f"+{int(row['marginal'])} from {row['design_id']}"
            )


@main.command("export-research")
@_IN
@_CONTEXTS
@click.option(
    "--clinvar-release", required=True, help="ClinVar release id for provenance, e.g. 20260715."
)
@click.option(
    "--commit", default="", help="Commit SHA for provenance; read from git HEAD if omitted."
)
@click.option(
    "--aenmd",
    "aenmd_table",
    type=click.Path(exists=True, dir_okay=False, path_type=Path),
    help="The aenmd verdicts table, to add the model tier's agreement with the rule tier.",
)
@click.option(
    "--nmdetective",
    "nmdetective_table",
    type=click.Path(exists=True, dir_okay=False, path_type=Path),
    help="The NMDetective-AI table, to add the deep model's efficiency separation to the atlas.",
)
@_OUT
def export_research(
    table: Path,
    contexts: Path,
    clinvar_release: str,
    commit: str,
    aenmd_table: Path | None,
    nmdetective_table: Path | None,
    out: Path,
) -> None:
    """Build the researcher dashboard aggregate: coverage frontiers and per-disease coverage.

    Reads the normalized disease table and the contexts table; writes one small JSON of aggregates,
    never per-variant rows, so the dashboard ships without loading the variant set.
    """

    if not commit:
        try:
            commit = subprocess.run(
                ["git", "rev-parse", "--short", "HEAD"],
                capture_output=True,
                text=True,
                check=True,
            ).stdout.strip()
        except (subprocess.CalledProcessError, FileNotFoundError):
            commit = ""

    aggregate = build_research_aggregate(
        read_table(table),
        read_table(contexts),
        clinvar_release=clinvar_release,
        commit=commit,
        aenmd=read_table(aenmd_table) if aenmd_table is not None else None,
        nmdetective=read_table(nmdetective_table) if nmdetective_table is not None else None,
    )
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(aggregate.to_json())
    click.echo(
        f"research aggregate: {aggregate.provenance['condition_entities']} condition entities over "
        f"{aggregate.provenance['qualifying_variants']} variants, "
        f"{len(aggregate.condition_coverage_top)} in the top list"
    )


@main.command("contexts")
@click.argument("variants", type=click.Path(exists=True, dir_okay=False, path_type=Path))
@click.option("--annotation", required=True, type=click.Path(exists=True, path_type=Path))
@click.option("--transcripts", required=True, type=click.Path(exists=True, path_type=Path))
@click.option("--proteins", required=True, type=click.Path(exists=True, path_type=Path))
@_OUT
def variant_contexts(
    variants: Path, annotation: Path, transcripts: Path, proteins: Path, out: Path
) -> None:
    """Read the stop-codon context of every variant on its gene's MANE Select transcript."""

    models = load_transcripts(annotation, transcripts)
    by_gene = {model.gene_id: model for model in models.values()}
    table = read_table(variants)

    disagreeing = disagreements_with_protein(table, by_gene, read_sequences(proteins))
    if disagreeing:
        listed = "\n".join(
            f"  {d.transcript_id} residue {d.protein_position}: "
            f"{d.codon} translates to {d.translated}, protein carries {d.in_protein}"
            for d in disagreeing[:10]
        )
        raise click.ClickException(
            f"{len(disagreeing)} placements disagree with the reference protein\n{listed}"
        )

    found = contexts_for(table, by_gene)
    out.parent.mkdir(parents=True, exist_ok=True)
    write_table(found, out)
    scoreable = int(found["scoreable"].sum())
    click.echo(f"{scoreable} of {len(found)} variants have a scoreable context")
    for reason, count in found.loc[~found["scoreable"], "reason"].value_counts().items():
        click.echo(f"  {count} {reason}")


@main.command("score")
@click.argument("contexts", type=click.Path(exists=True, dir_okay=False, path_type=Path))
@click.argument(
    "training",
    nargs=-1,
    required=True,
    type=click.Path(exists=True, dir_okay=False, path_type=Path),
)
@_OUT
def score_contexts(contexts: Path, training: tuple[Path, ...], out: Path) -> None:
    """Predict readthrough for every context in CONTEXTS, under each measured therapy.

    A therapy is named after its TRAINING table, so `features_G418.tsv.gz` scores as G418.
    """

    variants = read_table(contexts)
    scoreable = variants[variants["scoreable"]].set_index("variant_id")
    scored = []

    # Sorted so the output does not depend on the order the caller's shell or workflow globbed them.
    for path in sorted(training):
        therapy = path.name.split(".")[0].removeprefix("features_")
        model = fit(read_table(path))
        known = model.known_levels(scoreable)
        predicted = model.predict(scoreable[known])

        held_out = path.with_name(f"predictions_{therapy}.tsv.gz")
        if not held_out.exists():
            raise click.ClickException(
                f"{held_out} is absent; a prediction without its held-out error has no uncertainty"
            )
        spread = relative_error_quantile(read_table(held_out))

        rows = scoreable.assign(
            therapy_id=therapy,
            readthrough_predicted=predicted,
            readthrough_low=(predicted * (1 - spread)).clip(lower=0.0),
            readthrough_high=predicted * (1 + spread),
        )
        rows["status"] = pd.Series(
            ["present" if seen else "missing" for seen in known], index=scoreable.index
        )
        rows["reason"] = ["" if seen else "not_available" for seen in known]
        scored.append(rows.reset_index())
        click.echo(
            f"{therapy:<12} scored {int(known.sum())} of {len(scoreable)} contexts; "
            f"95% of held-out predictions landed within {spread * 100:.0f}% of the value"
        )

    table = pd.concat(scored, ignore_index=True)
    out.parent.mkdir(parents=True, exist_ok=True)
    write_table(
        table[
            [
                "variant_id",
                "gene_symbol",
                "transcript_id",
                "protein_position",
                "stop_type",
                "therapy_id",
                "readthrough_predicted",
                "readthrough_low",
                "readthrough_high",
                "status",
                "reason",
                "review_stars",
            ]
        ],
        out,
    )
    click.echo(f"wrote {len(table)} variant by therapy rows to {out}")


@main.command("nmd")
@_IN
@_OUT
def nmd(table: Path, out: Path) -> None:
    """Two rule-based NMD escape predictors per variant, and where they disagree.

    The guideline predictor is the 50-nt last-junction rule (the one ClinGen PVS1 uses); the full
    rule set adds the start-proximal and long-exon escapes. Disagreement is where the fuller rules
    escape a stop the guideline calls decay — not evidence of which rule is right (ADR-0016).
    """

    predictors = nmd_predictors(read_table(table))
    out.parent.mkdir(parents=True, exist_ok=True)
    write_table(predictors, out)

    atlas = disagreement_atlas(predictors)
    click.echo(f"{atlas['scoreable']} scoreable stops")
    click.echo(
        f"  guideline (50-nt) escape: {atlas['escape_guideline']:>6} "
        f"({atlas['guideline_fraction'] * 100:>4.1f}%)"
    )
    click.echo(
        f"  full rule set escape:     {atlas['escape_full_rules']:>6} "
        f"({atlas['full_rules_fraction'] * 100:>4.1f}%)"
    )
    click.echo(
        f"  predictors disagree:      {atlas['disagree']:>6} "
        f"({atlas['disagree_fraction'] * 100:>4.1f}%)"
    )
    click.echo(
        f"    driven by start-proximal {atlas['driven_by_start_proximal']}, "
        f"long-exon {atlas['driven_by_long_exon']}, both {atlas['driven_by_both']}"
    )


@main.command("aenmd-verdicts")
@click.option(
    "--aenmd",
    "aenmd_rules",
    required=True,
    type=click.Path(exists=True, dir_okay=False, path_type=Path),
    help="aenmd's per-transcript rule table (results/aenmd.tsv), from scripts/aenmd_nmd.R.",
)
@click.option(
    "--nonsense",
    required=True,
    type=click.Path(exists=True, dir_okay=False, path_type=Path),
    help="The ClinVar nonsense variants, for the coordinates the aenmd key is built from.",
)
@click.option(
    "--mane-summary",
    required=True,
    type=click.Path(exists=True, dir_okay=False, path_type=Path),
    help="The MANE summary, for the one-to-one RefSeq/Ensembl transcript pairing.",
)
@click.option(
    "--contexts",
    "contexts_table",
    type=click.Path(exists=True, dir_okay=False, path_type=Path),
    help="The context table; when given, reports how the rule tier and aenmd agree.",
)
@_OUT
def aenmd_verdicts_command(
    aenmd_rules: Path, nonsense: Path, mane_summary: Path, contexts_table: Path | None, out: Path
) -> None:
    """aenmd's NMD escape verdict per variant, on its MANE Select transcript (ADR-0017).

    aenmd is the model tier's real, independent rule implementation. This reads its per-transcript
    output down to one verdict per variant on the MANE transcript, keeping every variant aenmd did
    not score marked unavailable with the reason. With `--contexts`, it reports agreement with the
    rule tier — the check that the hand-rolled `full_rules` reproduces the published tool.
    """

    rules = read_aenmd_rules(str(aenmd_rules))
    mane = mane_ensembl_by_gene(str(mane_summary))
    nons = pd.read_csv(nonsense, sep="\t", dtype={"chrom": str})
    verdicts = aenmd_verdicts(rules, nons, mane)
    out.parent.mkdir(parents=True, exist_ok=True)
    write_table(verdicts, out)

    available = int(verdicts["aenmd_available"].sum())
    click.echo(
        f"{len(verdicts):,} variants, aenmd available for {available:,} "
        f"({available / len(verdicts) * 100:.1f}%) → {out}"
    )
    if contexts_table is not None:
        agree = model_agreement(nmd_predictors(read_table(contexts_table)), verdicts)
        if agree.get("both_available"):
            click.echo(
                f"  full_rules vs aenmd agree: {agree['full_rules_vs_aenmd_agree']:,}/"
                f"{agree['both_available']:,} "
                f"({agree['full_rules_vs_aenmd_agree_fraction'] * 100:.2f}%); "
                f"aenmd-only escapes {agree['full_decay_aenmd_escape']}, "
                f"rule-only escapes {agree['full_escape_aenmd_decay']}"
            )


@main.command("trna-coverage")
@click.argument("contexts", type=click.Path(exists=True, dir_okay=False, path_type=Path))
@_OUT
def trna_coverage(contexts: Path, out: Path) -> None:
    """Rank suppressor tRNA designs by the pathogenic variants each would reach."""

    table = read_table(contexts)
    scoreable = table[table["scoreable"]]
    coverage = coverage_by_design(scoreable)
    out.parent.mkdir(parents=True, exist_ok=True)
    write_table(coverage, out)

    click.echo(f"{len(coverage)} designs over {len(scoreable)} variants")
    for row in coverage.head(5).itertuples():
        click.echo(
            f"  {row.design_id:<7} {row.conservative:>6} conservative, "
            f"{row.restores_exactly:>6} restored exactly"
        )


@main.command("trna-panel")
@click.argument("contexts", type=click.Path(exists=True, dir_okay=False, path_type=Path))
@click.option(
    "--objective",
    type=click.Choice(["variants", "genes"]),
    default="variants",
    help="What each added design should cover the most of.",
)
@_OUT
def trna_panel(contexts: Path, objective: str, out: Path) -> None:
    """Greedy suppressor-tRNA panel: the fewest designs covering the most variants or genes.

    Coverage is exact restoration — the design puts the native residue back at the stop — and
    nothing more; it is model coverage, not a therapeutic or clinical claim, with no safety axis.
    """

    table = read_table(contexts)
    frontier = coverage_frontier(table, objective=objective)
    out.parent.mkdir(parents=True, exist_ok=True)
    write_table(frontier, out)

    total = int(frontier["cumulative"].iloc[-1]) if len(frontier) else 0
    click.echo(
        f"{len(frontier)} exact-restoration designs observed, covering all {total} {objective} "
        f"in the qualifying ClinVar nonsense-variant set"
    )
    for k in (1, 3, 5, 10):
        if k <= len(frontier):
            row = frontier.iloc[k - 1]
            click.echo(
                f"  panel of {k:>2}: {int(row['cumulative']):>5} of {total} {objective} "
                f"({row['cumulative_fraction'] * 100:>5.1f}%), "
                f"+{int(row['marginal'])} from {row['design_id']}"
            )


@main.command("extensions")
@click.argument("transcripts", type=click.Path(exists=True, dir_okay=False, path_type=Path))
@click.option(
    "--annotation",
    required=True,
    type=click.Path(exists=True, dir_okay=False, path_type=Path),
    help="Transcript lengths, as the calibration step writes them.",
)
@_OUT
def extensions(transcripts: Path, annotation: Path, out: Path) -> None:
    """Find how far past each native stop the next in-frame stop lies."""

    windows = extension_windows(transcripts, read_table(annotation))
    # Written aside and moved into place, so a partial file is never mistaken for a finished one.
    staging = out.with_suffix(out.suffix + ".tmp")
    write_table(windows, staging)
    staging.replace(out)
    usable = int(windows["extension"].notna().sum())
    click.echo(f"{len(windows):,} coding transcripts, {usable:,} with a next in-frame stop")


@main.command("readthrough")
@click.argument("counts", type=click.Path(exists=True, dir_okay=False, path_type=Path))
@click.option("--gtf", required=True, type=click.Path(exists=True, path_type=Path))
@click.option(
    "--samplesheet",
    required=True,
    type=click.Path(exists=True, dir_okay=False, path_type=Path),
    help="The staged runs, naming each library's treatment and replicate.",
)
@click.option("--treated", required=True, help="The treatment arm being tested.")
@click.option("--control", required=True, help="The arm it is tested against.")
@click.option("--dataset", required=True, help="The experiment whose libraries form the contrast.")
@click.option(
    "--manifest",
    required=True,
    type=click.Path(exists=True, dir_okay=False, path_type=Path),
    help="The calibration manifest for this dataset. Refused when it records a failure.",
)
@click.option(
    "--published-lengths",
    type=(int, int),
    help="Run the sensitivity arm on this inclusive length range instead of the selected set. "
    "The manifest is still required and must still pass.",
)
@click.option(
    "--paired",
    is_flag=True,
    help="Compare within replicate. Only where the metadata documents which libraries were "
    "prepared together; matching replicate numbers across arms do not.",
)
@_OUT
def readthrough_assay(
    counts: Path,
    gtf: Path,
    samplesheet: Path,
    treated: str,
    control: str,
    dataset: str,
    manifest: Path,
    published_lengths: tuple[int, int] | None,
    paired: bool,
    out: Path,
) -> None:
    """Test COUNTS for the readthrough signature, one contrast at a time."""

    # Nothing else asks whether these libraries are ribosome profiles. A dataset that failed its
    # predeclared calibration has no result to report, and the way that goes wrong is a contrast
    # computed quietly on libraries nobody looked at.
    try:
        calibration = read_manifest(manifest)
    except (ValueError, KeyError, OSError) as error:
        raise click.ClickException(str(error)) from error
    if calibration.dataset != dataset:
        raise click.ClickException(f"{manifest} calibrates {calibration.dataset}, not {dataset}")
    # The sensitivity arm still passes through the calibration gate; only which lengths it sums
    # differs, so it can never be run on a dataset whose libraries failed.
    lengths_used = (
        tuple(range(published_lengths[0], published_lengths[1] + 1))
        if published_lengths is not None
        else calibration.lengths
    )

    runs = _validated(StagedRuns, read_table(samplesheet), samplesheet)
    arms = runs.loc[
        (runs["assay"] == "riboseq")
        & (runs["dataset"] == dataset)
        & (runs["treatment"].isin([treated, control])),
        ["sample", "treatment", "replicate"],
    ]
    if arms.empty:
        raise click.ClickException(f"{samplesheet} names no {dataset} libraries for that contrast")

    measured = read_table(counts)
    # Counts arrive stratified by footprint length so that one pass over the alignments serves both
    # the selected set and the published window. Which lengths this contrast uses comes from the
    # manifest, never from an argument, so it cannot differ from what the manifest records.
    if "length" in measured.columns:
        measured = collapse_lengths(measured, lengths_used)
    # A library named for the contrast but absent from the counts would quietly shrink the
    # comparison — three against three becoming two against three without a word.
    if absent := sorted(set(arms["sample"]) - set(measured["sample"])):
        raise click.ClickException(f"{counts} has no rows for: {', '.join(absent)}")
    kept = qualifying(
        measured,
        transcript_genes(gtf),
        overlapping_downstream_cds(gtf),
        samples=set(arms["sample"]),
    )
    lengths = ", ".join(str(length) for length in lengths_used)
    click.echo(f"{len(measured):,} rows, {kept['transcript'].nunique():,} transcripts qualifying")
    arm = "published window" if published_lengths is not None else "selected set"
    click.echo(f"{arm}: {lengths} nt")

    ratios = library_ratios(kept)
    compare = paired_effect if paired else unpaired_effect
    try:
        effects = {
            quantity: compare(ratios, quantity, arms, treated, control)
            for quantity in ("downstream_occupancy", "termination_occupancy", "frame_gap")
        }
    except ValueError as error:
        raise click.ClickException(str(error)) from error

    met = signature(effects)
    summary = pd.DataFrame(
        {
            "quantity": list(effects),
            "mean_difference": [e.mean_difference for e in effects.values()],
            "ci_low": [e.interval[0] for e in effects.values()],
            "ci_high": [e.interval[1] for e in effects.values()],
            "consistent": [e.consistent for e in effects.values()],
        }
    )
    write_table(summary, out)
    # The per-library figures the effects were built from, so a pooled proportion resting on a
    # handful of transcripts is visible rather than buried under three summary rows.
    diagnostics = out.with_name(f"{out.stem}_by_library{out.suffix}")
    write_table(ratios.merge(arms, on="sample"), diagnostics)

    click.echo(f"{treated} against {control}, {'paired' if paired else 'unpaired'}")
    for row in summary.itertuples():
        click.echo(
            f"  {row.quantity:>22}: {row.mean_difference:+.4g} "
            f"[{row.ci_low:+.4g}, {row.ci_high:+.4g}] "
            f"{'consistent' if row.consistent else 'inconsistent'}"
        )
    for condition, held in met.items():
        click.echo(f"  {condition:>22}: {'yes' if held else 'no'}")
    click.echo(f"  {'signature':>22}: {'complete' if all(met.values()) else 'incomplete'}")
    if stalling(effects):
        click.echo(f"  {'stalling':>22}: yes — raised at the stop, not beyond it")


@main.command("landscape")
@click.argument("contexts", type=click.Path(exists=True, dir_okay=False, path_type=Path))
@click.argument("scores", type=click.Path(exists=True, dir_okay=False, path_type=Path))
@_OUT
@click.option(
    "--summary", type=click.Path(dir_okay=False, path_type=Path), help="Where to write the counts."
)
def amenability_landscape(contexts: Path, scores: Path, out: Path, summary: Path | None) -> None:
    """Bring transcript survival, readthrough and residue tolerance together, variant by variant."""

    table = landscape(read_table(contexts), read_table(scores))
    out.parent.mkdir(parents=True, exist_ok=True)
    write_table(table, out)

    counts = summarise(table, Thresholds())
    if summary is not None:
        write_table(counts, summary)

    click.echo(f"{len(table)} variants placed")
    tolerable = int((table["tolerable_insertion_share"] >= TOLERABLE_SHARE).sum())
    click.echo(
        f"  {int(table['escapes_decay_by_rule'].sum())} expected to escape decay; "
        f"{tolerable} where most insertions are tolerable"
    )
    for row in counts.to_dict("records"):
        click.echo(
            f"  readthrough >= {row['readthrough_threshold'] * 100:.1f}%: "
            f"{row['reaches_threshold']:>6} reach it, "
            f"{row['all_conditions']:>5} meet every condition "
            f"({row['all_conditions_lower_bound']} on the interval's lower bound)"
        )


def _validated(
    schema: type[pandera.pandas.DataFrameModel], frame: pd.DataFrame, source: Path
) -> pd.DataFrame:
    try:
        return schema.validate(frame, lazy=True)
    except pandera.errors.SchemaErrors as error:
        raise click.ClickException(f"{source} does not match {schema.__name__}\n{error}") from error
