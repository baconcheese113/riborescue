#!/usr/bin/env nextflow

/*
 * RiboRescue — three workflows, selected with --step.
 *
 *   --step amenability   ClinVar's pathogenic nonsense variants, placed on their MANE Select
 *                        transcripts, scored under each therapy, ranked by suppressor design, and
 *                        joined into the landscape of what is plausibly addressable
 *   --step train         fits the readthrough model from measured labels
 *   --step score         triages submitted variants against the upstream handoff
 *
 * Only the workflows consuming Ribo-seq output require the upstream handoff; the amenability path
 * runs from public annotation alone. Every step calls the installed riborescue command; no analysis
 * logic lives in Groovy.
 */

include { AMENABILITY_LANDSCAPE } from './modules/local/amenability_landscape.nf'
include { CLINVAR_VARIANTS } from './modules/local/clinvar_variants.nf'
include { SCORE_VARIANTS   } from './modules/local/score_variants.nf'
include { TRIAGE_VARIANTS  } from './modules/local/triage_variants.nf'
include { TRNA_COVERAGE    } from './modules/local/trna_coverage.nf'
include { VALIDATE_HANDOFF } from './modules/local/validate_handoff.nf'
include { VALIDATE_LABELS  } from './modules/local/validate_labels.nf'
include { VARIANT_CONTEXTS } from './modules/local/variant_contexts.nf'

workflow HANDOFF {
    main:
    if( !params.handoff || !params.results_root )
        error "--handoff and --results_root are required: the manifest naming the " +
              "${params.upstream_pipeline} outputs to consume, and the tree holding them"

    VALIDATE_HANDOFF(
        channel.fromPath(params.handoff, checkIfExists: true),
        file(params.results_root, checkIfExists: true)
    )
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

workflow TRAIN {
    main:
    if( !params.labels )
        error '--labels is required: the measured readthrough efficiency table to fit against'

    HANDOFF()
    VALIDATE_LABELS(channel.fromPath(params.labels, checkIfExists: true))
}

workflow SCORE {
    main:
    if( !params.variants )
        error '--variants is required: the variant table to score'

    HANDOFF()
    TRIAGE_VARIANTS(channel.fromPath(params.variants, checkIfExists: true))
}

workflow {
    main:
    if( params.step == 'amenability' )
        AMENABILITY()
    else if( params.step == 'train' )
        TRAIN()
    else if( params.step == 'score' )
        SCORE()
    else
        error "--step must be 'amenability', 'train' or 'score', not '${params.step}'"
}
