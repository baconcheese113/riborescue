// The shape of riborescue.json, written by `riborescue export-web`. Kept beside the app so a schema
// change surfaces here as a type error rather than as a blank cell.

export interface Therapy {
  id: string;
  // The compound the id names. "SRI" is SRI-41315, which is not the SRI-37240 of the safety
  // control's Ribo-seq dataset, so the fuller name is what the interface shows.
  name: string;
  available: boolean;
  readthrough: number | null;
  low: number | null;
  high: number | null;
  reason: string | null;
}

export interface Suppressor {
  design: string;
  restores_exactly: boolean;
}

// Geometric base-editing reachability of one variant. A candidate placement, never eligibility:
// not editing efficiency, off-target activity, delivery, tissue access, or splice.
export interface Editing {
  reach_class:
    | "base_editable_exact"
    | "base_editable_alternative"
    | "not_base_editable_under_panel";
  reachable: boolean;
  editor: string | null;
  strand: string | null;
  restores: string | null;
  window_position: number | null;
  bystander_free: boolean | null;
  candidate_guides: number;
}

// The escape-map denominator flow: every variant accounted for, none dropped silently. Percentages
// are computed against `scoreable`; the total and unscoreable count travel with it.
export interface EscapeSummary {
  panel: string;
  total: number;
  scoreable: number;
  unscoreable: number;
  exact: number;
  alternative: number;
  not_editable: number;
  reachable: number;
  reachable_bystander_free: number;
  exact_bystander_free: number;
  caveat: string;
}

export interface Variant {
  id: string;
  gene: string;
  protein_position: number | null;
  stop: string;
  residue: string;
  review_stars: number | null;
  escapes_decay: boolean;
  tolerable_insertion_share: number;
  best: { therapy: string; readthrough: number; low: number } | null;
  therapies: Therapy[];
  suppressor: Suppressor | null;
  nmd: NmdVerdict | null;
  editing: Editing | null;
}

export interface AenmdVerdict {
  available: boolean;
  escape: boolean;
  reason: string;
}

export interface NmdetectiveVerdict {
  available: boolean;
  efficiency: number | null;
  reason: string;
}

export interface NmdVerdict {
  escape_guideline: boolean;
  escape_full_rules: boolean;
  disagree: boolean;
  rules: {
    last_exon: boolean;
    within_last_junction: boolean;
    start_proximal: boolean;
    long_exon: boolean;
  };
  aenmd?: AenmdVerdict | null;
  nmdetective?: NmdetectiveVerdict | null;
}

export interface SafetyAtlas {
  measured_therapy: string;
  cell_line: string;
  canonical_stops_scored: number;
  analysed: number;
  unavailable_reason: string;
  concordance: { rho: number; low: number; high: number };
  concordance_label: string;
  quadrants: Record<string, number>;
  per_therapy: Record<string, string>;
  caveat: string;
}

export interface WebTable {
  contracts_version: string;
  status: {
    readthrough_control: string;
    readthrough_detail: string;
    safety_atlas: string;
    safety_detail: string;
  };
  therapies: string[];
  therapy_names: Record<string, string>;
  safety: SafetyAtlas | null;
  escape: EscapeSummary | null;
  variants: Variant[];
}

// The shape of riborescue_research.json, written by `riborescue export-research`. Aggregates only —
// coverage frontiers and per-disease coverage — never the variant rows.

export interface FrontierStep {
  rank: number;
  design_id: string;
  cumulative: number;
  cumulative_fraction: number;
  marginal: number;
}

export interface ConditionCoverage {
  condition_name: string;
  medgen: string;
  omim: string;
  orphanet: string;
  eligible_variants: number;
  model_covered: number;
  covered_fraction: number;
  reach: boolean;
  complete: boolean;
  genes: number;
  designs_contributing: number;
  mapping_completeness: string;
}

// One proposed programme. `direct_scope` is what the experiment actually puts under measurement;
// `potential_*` is who would benefit if the result generalises, and the two are never merged.
export interface Experiment {
  experiment_id: string;
  question: string;
  why_it_matters: string;
  what_the_lab_does: string;
  comparison: string;
  success_criterion: string;
  if_it_fails: string;
  evidence_gap_reason: string;
  direct_scope: string;
  generalisation_required: string;
  assay: string;
  model_system: string;
  endpoint: string;
  decision_rule: string;
  replicates: string;
  replicate_endpoint: string;
  replicate_effect: string;
  replicate_variance_source: string;
  replicate_design: string;
  replicate_alpha_power: string;
  replicate_method: string;
  evidence_grade: string;
  complexity: string;
  safety_relevant: boolean;
  provenance: string;
  claims_named: string;
  resolves: string;
  reach_rule: string;
  evidence_gap: number;
  feasibility: number;
  potential_variants: number;
  potential_genes: number;
  potential_conditions: number;
  open_questions_addressed: number;
  on_frontier: boolean;
  dominated_by: string;
}

export interface ResearchAggregate {
  provenance: {
    clinvar_release: string;
    commit: string;
    qualifying_variants: number;
    scoreable_variants: number;
    condition_entities: number;
  };
  mapping_completeness: Record<string, number>;
  escape?: EscapeSummary | null;
  reach_denominator: {
    eligible_condition_entities: number;
    unique_variant_condition_pairs: number;
    excluded_placeholder_or_unmapped_rows: number;
    reachable_entities: number;
    note: string;
  };
  nmd: {
    scoreable: number;
    escape_guideline: number;
    escape_full_rules: number;
    guideline_fraction: number;
    full_rules_fraction: number;
    disagree: number;
    disagree_fraction: number;
    driven_by_start_proximal: number;
    driven_by_long_exon: number;
    driven_by_both: number;
    aenmd?: {
      scoreable: number;
      aenmd_available: number;
      both_available: number;
      aenmd_available_fraction: number;
      aenmd_escape: number;
      aenmd_escape_fraction: number;
      full_rules_vs_aenmd_agree: number;
      full_rules_vs_aenmd_agree_fraction: number;
      full_escape_aenmd_decay: number;
      full_decay_aenmd_escape: number;
      guideline_vs_aenmd_agree_fraction: number;
    };
    nmdetective?: {
      scoreable: number;
      available: number;
      both_available: number;
      available_fraction: number;
      efficiency_median?: number;
      guideline?: { escape_mean_efficiency: number | null; decay_mean_efficiency: number | null; separation: number | null };
      full_rules?: { escape_mean_efficiency: number | null; decay_mean_efficiency: number | null; separation: number | null };
    };
  };
  frontiers: {
    variants: FrontierStep[];
    genes: FrontierStep[];
    conditions: FrontierStep[];
  };
  condition_coverage_top: ConditionCoverage[];
  experiments: Experiment[];
  caveats: Record<string, string>;
}

export interface LandscapeThreshold {
  readthrough_threshold: number;
  variants: number;
  escapes_decay: number;
  reaches_threshold: number;
  reaches_threshold_lower_bound: number;
  most_insertions_tolerable: number;
  all_conditions: number;
  all_conditions_lower_bound: number;
}


/** The evidence payload: the controls, the calibration, and the checks behind them. */
export interface EvidenceQuantity {
  quantity: string;
  mean_difference: number | null;
  ci_low: number | null;
  ci_high: number | null;
  consistent: boolean;
}

export interface EvidenceLibrary {
  sample: string;
  treatment: string | null;
  transcripts: number;
  downstream_occupancy?: number | null;
  termination_occupancy?: number | null;
  frame_gap?: number | null;
}

export interface EvidenceContrast {
  quantities: EvidenceQuantity[];
  libraries: EvidenceLibrary[];
}

export interface FrameByLength {
  length: number;
  frame0: number;
  frame1: number;
  frame2: number;
  frame0_share: number | null;
  library_share: number | null;
  kept: boolean;
}

export interface EvidenceCalibration {
  dataset: string;
  lengths: number[];
  surveyed: number[];
  passes: boolean;
  libraries: {
    sample: string;
    psites: number;
    frame0_share: number;
    dominant_length: number;
    offset_from_5: number;
    failures: string[];
  }[];
  frame_by_length?: FrameByLength[];
}

export interface PeriodicityPoint {
  region: "start" | "stop";
  treatment: string | null;
  distance: number;
  scaled: number | null;
  libraries: number;
}

export interface CodonOccupancy {
  codon: string;
  amino_acid: string;
  site: "a" | "p";
  occupancy: number | null;
  occupancy_sd: number | null;
  libraries: number;
}

export interface KineticsNull {
  permutations_completed: number;
  permutations_required: number;
  analysis_status: "complete" | "incomplete";
  resolution: number | null;
  rows: KineticsNullRow[];
}

export interface KineticsNullRow {
  drug: string;
  shuffle: string;
  gain: number | null;
  null_mean: number | null;
  null_sd: number | null;
  null_max: number | null;
  p_familywise: number | null;
  permutations: number;
}

export interface SafetyConcordance {
  rho: number;
  low: number;
  high: number;
  analysed: number;
  canonical_stops_scored: number;
  quadrants: Record<string, number>;
  points: { gene: string; predicted: number | null; measured: number | null; group: string }[];
}

export interface ModelParity {
  drug: string;
  r2_mean: number | null;
  r2_sd: number | null;
  rounds: number;
  ceiling: number | null;
}

export interface EvidenceTable {
  contracts_version: string;
  provenance: { dataset: string; commit: string; scope: string };
  readthrough: Record<string, EvidenceContrast> | null;
  calibration: EvidenceCalibration | null;
  periodicity: PeriodicityPoint[] | null;
  codon_occupancy: CodonOccupancy[] | null;
  kinetics_null: KineticsNull | null;
  safety: SafetyConcordance | null;
  model_parity: ModelParity[] | null;
}
