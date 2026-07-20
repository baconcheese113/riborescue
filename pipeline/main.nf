#!/usr/bin/env nextflow

/*
 * RiboRescue — two workflows over one validated upstream handoff.
 *
 *   --step train   fits the readthrough model from measured labels
 *   --step score   applies it to submitted variants
 *
 * Every step calls the installed riborescue command; no analysis logic lives in Groovy.
 */

include { VALIDATE_HANDOFF } from './modules/local/validate_handoff.nf'
include { VALIDATE_LABELS  } from './modules/local/validate_labels.nf'
include { TRIAGE_VARIANTS  } from './modules/local/triage_variants.nf'

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
    if( params.step == 'train' )
        TRAIN()
    else if( params.step == 'score' )
        SCORE()
    else
        error "--step must be 'train' or 'score', not '${params.step}'"
}
