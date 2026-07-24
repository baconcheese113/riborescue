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

export interface DiseaseCoverage {
  disease_name: string;
  medgen: string;
  omim: string;
  orphanet: string;
  eligible_variants: number;
  model_covered: number;
  covered_fraction: number;
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
    diseases: number;
  };
  mapping_completeness: Record<string, number>;
  frontiers: {
    variants: FrontierStep[];
    genes: FrontierStep[];
    diseases: FrontierStep[];
  };
  disease_coverage_top: DiseaseCoverage[];
  caveats: Record<string, string>;
}
