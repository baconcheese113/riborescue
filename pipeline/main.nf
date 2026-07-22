#!/usr/bin/env nextflow

/*
 * RiboRescue — three workflows, selected with --step.
 *
 *   --step amenability   ClinVar's pathogenic nonsense variants, placed on their MANE Select
 *                        transcripts, scored under each therapy, ranked by suppressor design, and
 *                        joined into the landscape of what is plausibly addressable
 *   --step reads         quality control, adapter trimming and quality control again over the
 *                        staged Ribo-seq and RNA-seq runs
 *   --step train         fits the readthrough model from measured labels
 *   --step score         triages submitted variants against the upstream handoff
 *
 * Only the workflows consuming Ribo-seq output require the upstream handoff; the amenability path
 * runs from public annotation alone. Every step calls the installed riborescue command; no analysis
 * logic lives in Groovy.
 */

include { ALIGNMENT_SUMMARY } from './modules/local/alignment_summary.nf'
include { AMENABILITY_LANDSCAPE } from './modules/local/amenability_landscape.nf'
include { CLINVAR_VARIANTS } from './modules/local/clinvar_variants.nf'
include { CONTAMINANT_INDEX    } from './modules/local/contaminant_index.nf'
include { DEPLETE_CONTAMINANTS } from './modules/local/deplete_contaminants.nf'
include { STAR_ALIGN       } from './modules/local/star_align.nf'
include { STAR_INDEX       } from './modules/local/star_index.nf'
include { CUTADAPT         } from './modules/local/cutadapt.nf'
include { FASTQC as FASTQC_RAW     } from './modules/local/fastqc.nf'
include { FASTQC as FASTQC_TRIMMED } from './modules/local/fastqc.nf'
include { MULTIQC          } from './modules/local/multiqc.nf'
include { TRIM_SUMMARY     } from './modules/local/trim_summary.nf'
include { SCORE_VARIANTS   } from './modules/local/score_variants.nf'
include { TRIAGE_VARIANTS  } from './modules/local/triage_variants.nf'
include { TRNA_COVERAGE    } from './modules/local/trna_coverage.nf'
include { VALIDATE_HANDOFF } from './modules/local/validate_handoff.nf'
include { VALIDATE_LABELS  } from './modules/local/validate_labels.nf'
include { VARIANT_CONTEXTS } from './modules/local/variant_contexts.nf'

// A read named relatively is found beside the sheet naming it, so a staged sheet stays valid
// wherever it is mounted.
def locate(sheet, String read) {
    read.startsWith('/') ? file(read) : sheet.parent.resolve(read)
}

workflow HANDOFF {
    main:
    if( !params.handoff || !params.results_root )
        error "--handoff and --results_root are required: the manifest naming the " +
              "${params.upstream_pipeline} outputs to consume, and the tree holding them"

    VALIDATE_HANDOFF(
        channel.fromPath(params.handoff, checkIfExists: true),
        file(params.results_root, checkIfExists: true)
    )

    emit:
    VALIDATE_HANDOFF.out.manifest
}

workflow AMENABILITY {
    main:
    def required = [
        clinvar: params.clinvar,
        mane_annotation: params.mane_annotation,
        mane_transcripts: params.mane_transcripts,
        mane_proteins: params.mane_proteins,
        training: params.training,
        held_out: params.held_out,
    ]
    def missing = required.findAll { _name, value -> !value }.keySet()
    if( missing )
        error "missing required inputs: ${missing.join(', ')} — fetch them with `riborescue fetch`"

    CLINVAR_VARIANTS(channel.fromPath(params.clinvar, checkIfExists: true))
    VARIANT_CONTEXTS(
        CLINVAR_VARIANTS.out.variants,
        file(params.mane_annotation, checkIfExists: true),
        file(params.mane_transcripts, checkIfExists: true),
        file(params.mane_proteins, checkIfExists: true)
    )
    SCORE_VARIANTS(
        VARIANT_CONTEXTS.out.contexts,
        channel.fromPath(params.training, checkIfExists: true).collect(),
        channel.fromPath(params.held_out, checkIfExists: true).collect()
    )
    TRNA_COVERAGE(VARIANT_CONTEXTS.out.contexts)
    AMENABILITY_LANDSCAPE(VARIANT_CONTEXTS.out.contexts, SCORE_VARIANTS.out.scored)
}

workflow READS {
    main:
    if( !params.samplesheet )
        error '--samplesheet is required: the staged runs written by `riborescue stage-runs`'

    def sheet = file(params.samplesheet, checkIfExists: true)
    def runs = channel.of(sheet)
        .splitCsv(header: true, sep: '\t')
        .map { row ->
            def reads = row.layout == 'paired'
                ? [locate(sheet, row.fastq_1), locate(sheet, row.fastq_2)]
                : [locate(sheet, row.fastq_1)]
            tuple(row.subMap('sample', 'assay', 'layout', 'adapter_3p', 'adapter_3p_2',
                             'adapter_overlap', 'cut_5p'), reads)
        }

    FASTQC_RAW(runs, 'raw')
    CUTADAPT(runs)
    FASTQC_TRIMMED(CUTADAPT.out.reads, 'trimmed')
    TRIM_SUMMARY(CUTADAPT.out.report.collect(), sheet)

    // Alignment is skipped when no reference is named, so quality control and trimming stay
    // runnable without the hours an index costs.
    def aligned = params.genome && params.annotation
    if( aligned ) {
        // Structural RNA is most of a footprint library and aligns to many loci rather than none,
        // so it is removed before the genome sees the reads. A transcriptome library is a few
        // percent structural RNA and pays the whole depletion cost to remove almost nothing, so
        // only footprints are depleted.
        def to_align = CUTADAPT.out.reads
        if( params.contaminants && params.rdna ) {
            CONTAMINANT_INDEX(
                file(params.contaminants, checkIfExists: true),
                file(params.rdna, checkIfExists: true)
            )
            def by_assay = CUTADAPT.out.reads.branch { meta, _reads ->
                deplete: meta.assay == 'riboseq'
                keep: true
            }
            DEPLETE_CONTAMINANTS(by_assay.deplete, CONTAMINANT_INDEX.out.index)
            to_align = DEPLETE_CONTAMINANTS.out.reads.mix(by_assay.keep)
        }

        STAR_INDEX(
            file(params.genome, checkIfExists: true),
            file(params.annotation, checkIfExists: true)
        )
        STAR_ALIGN(to_align, STAR_INDEX.out.index)
        ALIGNMENT_SUMMARY(STAR_ALIGN.out.log.collect())
    }

    MULTIQC(
        FASTQC_RAW.out.zip
            .mix(
                FASTQC_TRIMMED.out.zip,
                CUTADAPT.out.report,
                aligned ? STAR_ALIGN.out.log : channel.empty()
            )
            .collect()
    )
}

workflow TRAIN {
    main:
    if( !params.labels )
        error '--labels is required: the measured readthrough efficiency table to fit against'

    // Combining with the validated manifest makes the handoff a prerequisite rather than a
    // neighbour: an invalid one stops the run before any label is read.
    HANDOFF()
    VALIDATE_LABELS(
        HANDOFF.out
            .combine(channel.fromPath(params.labels, checkIfExists: true))
            .map { _manifest, labels -> labels }
    )
}

workflow SCORE {
    main:
    if( !params.variants )
        error '--variants is required: the variant table to score'

    HANDOFF()
    TRIAGE_VARIANTS(
        HANDOFF.out
            .combine(channel.fromPath(params.variants, checkIfExists: true))
            .map { _manifest, variants -> variants }
    )
}

workflow {
    main:
    if( params.step == 'amenability' )
        AMENABILITY()
    else if( params.step == 'reads' )
        READS()
    else if( params.step == 'train' )
        TRAIN()
    else if( params.step == 'score' )
        SCORE()
    else
        error "--step must be 'amenability', 'reads', 'train' or 'score', not '${params.step}'"
}
