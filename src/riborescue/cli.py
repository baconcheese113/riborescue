"""The riborescue command line — every scientific step the pipeline runs enters through here."""

from pathlib import Path

import click
import pandas as pd
import pandera.errors
import pandera.pandas

from riborescue._version import __version__
from riborescue.contracts import Consequence
from riborescue.handoff import UpstreamHandoff
from riborescue.tables import ReadthroughLabels, TriageInput, TriageOutput, read_table, write_table
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


def _validated(
    schema: type[pandera.pandas.DataFrameModel], frame: pd.DataFrame, source: Path
) -> pd.DataFrame:
    try:
        return schema.validate(frame, lazy=True)
    except pandera.errors.SchemaErrors as error:
        raise click.ClickException(f"{source} does not match {schema.__name__}\n{error}") from error
