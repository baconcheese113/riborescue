"""The riborescue command line — every scientific step the pipeline runs enters through here."""

from pathlib import Path

import click
import pandas as pd
import pandera.errors
import pandera.pandas

from riborescue._version import __version__
from riborescue.baseline import fit, relative_error_quantile
from riborescue.clinvar import pathogenic_nonsense
from riborescue.context import contexts_for, disagreements_with_protein
from riborescue.contracts import Consequence, EvalConfig
from riborescue.evaluation import (
    SEED,
    UnsupportedEvalConfigError,
    bootstrap_ci,
    evaluate,
    split,
)
from riborescue.handoff import UpstreamHandoff
from riborescue.inputs import INPUTS, UnknownInputError, data_root, fetch
from riborescue.tables import (
    PathogenicNonsense,
    ReadthroughLabels,
    TriageInput,
    TriageOutput,
    read_table,
    write_table,
)
from riborescue.transcripts import load_transcripts, read_sequences
from riborescue.triage import classify, classify_table

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
    for path in features:
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

    for path in training:
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


def _validated(
    schema: type[pandera.pandas.DataFrameModel], frame: pd.DataFrame, source: Path
) -> pd.DataFrame:
    try:
        return schema.validate(frame, lazy=True)
    except pandera.errors.SchemaErrors as error:
        raise click.ClickException(f"{source} does not match {schema.__name__}\n{error}") from error
