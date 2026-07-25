// The shape of riborescue.json, written by `riborescue export-web`. Kept beside the app so a schema
// change surfaces here as a type error rather than as a blank cell.

export interface Therapy {
  id: string;
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
}

export interface AenmdVerdict {
  available: boolean;
  escape: boolean;
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
  status: {
    readthrough_control: string;
    readthrough_detail: string;
    safety_atlas: string;
    safety_detail: string;
  };
  therapies: string[];
  safety: SafetyAtlas | null;
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

export interface ResearchAggregate {
  provenance: {
    clinvar_release: string;
    commit: string;
    qualifying_variants: number;
    scoreable_variants: number;
    condition_entities: number;
  };
  mapping_completeness: Record<string, number>;
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
  };
  frontiers: {
    variants: FrontierStep[];
    genes: FrontierStep[];
    conditions: FrontierStep[];
  };
  condition_coverage_top: ConditionCoverage[];
  caveats: Record<string, string>;
}
