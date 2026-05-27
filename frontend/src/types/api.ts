export type Reference = {
  label: string
  url: string
  note: string
}

export type Motif = {
  id: string
  name: string
  factor: string
  source_type: "iupac" | "pwm"
  iupac_consensus: string | null
  pwm_matrix: Array<Record<"A" | "C" | "G" | "T", number>> | null
  organism_scope: string
  references: Reference[]
  limitations: string
}

export type BackgroundFrequencies = {
  A: number
  C: number
  G: number
  T: number
}

export type AnalysisSummary = {
  id: string
  motif_id: string
  motif_name: string
  sequence_count: number
  total_bases: number
  total_windows: number
  observed_hits: number
  expected_hits: number
  expected_hit_rate: number
  enrichment_test_count: number
  enrichment_model: "poisson_pwm_threshold"
  enrichment_p_value: number
  enrichment_q_value: number
  depletion_p_value: number
  two_sided_p_value: number
  fold_change: number
  background: BackgroundFrequencies
  status: "completed"
}

export type Hit = {
  id: number
  gene_id: string
  description: string
  strand: "+" | "-"
  start: number
  end: number
  matched_sequence: string
  score: number | null
  window_p_value: number
  gene_p_value: number
  q_value: number
}

export type AnalysisDetail = {
  summary: AnalysisSummary
  hits: Hit[]
}
