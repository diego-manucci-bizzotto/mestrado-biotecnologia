import {
  Activity,
  ArrowDownUp,
  CheckCircle2,
  ChevronRight,
  Dna,
  FileJson,
  FileSpreadsheet,
  Filter,
  Loader2,
  Microscope,
  Search,
  UploadCloud,
} from "lucide-react"
import { useMemo, useState } from "react"
import type { ChangeEvent, FormEvent } from "react"

import { Button } from "@/components/ui/button"

type Base = "A" | "C" | "G" | "T"
type RunMode = "pvalue" | "score"
type EngineMode = "custom" | "fimo" | "both"
type ViewTab = "overview" | "genes" | "matches" | "compare"
type SortMode = "pvalue" | "score" | "gene" | "position"
type StrandFilter = "all" | "forward" | "reverse"

type FastaProfile = {
  fileName: string
  sizeBytes: number
  records: number
  totalBases: number
  canonicalBases: number
  ambiguousBases: number
  gcRatio: number
  counts: Record<Base, number>
  firstHeader: string
}

type MotifMetadata = {
  filename: string
  motif_pattern: string
  motif_length: number
  score_threshold: number | null
  cutoff_mode: string
  background_frequencies: Record<Base, number>
}

type RawMotifMatch = {
  gene: string
  record_orientation: string
  strand: string
  strand_symbol: string
  source_strand: string
  source_strand_symbol: string
  matched_sequence: string
  position: number
  end_position: number
  source_position: number
  source_end_position: number
  score: number
  p_value: string
  p_value_numeric: number
  "-log10_pval": number | null
}

type ApiPayload = {
  metadata: MotifMetadata
  matches: RawMotifMatch[]
}

type ComparisonStats = {
  custom_total: number
  fimo_total: number
  overlap_total: number
  only_custom_total: number
  only_fimo_total: number
  overlap_ratio_custom: number
  overlap_ratio_fimo: number
}

type BothEnginesPayload = {
  engine: "both"
  custom: ApiPayload
  fimo: ApiPayload
  comparison: ComparisonStats
}

type SingleEnginePayload = {
  engine?: "custom" | "fimo"
  metadata: MotifMetadata
  matches: RawMotifMatch[]
}

type MotifMatch = RawMotifMatch & {
  negLog10PValue: number | null
}

type ScanResult = {
  engine: EngineMode
  metadata: MotifMetadata
  matches: MotifMatch[]
  custom?: {
    metadata: MotifMetadata
    matches: MotifMatch[]
  }
  fimo?: {
    metadata: MotifMetadata
    matches: MotifMatch[]
  }
  comparison?: ComparisonStats
  completedAt: string
}

type RunStats = {
  totalMatches: number
  uniqueGenes: number
  forwardSource: number
  reverseSource: number
  bestPValue: number | null
  maxScore: number | null
}

type GeneSummaryRow = {
  gene: string
  hits: number
  bestPValue: number
  bestScore: number
  firstPosition: number
  sourceStrands: string
}

const API_BASE = (import.meta.env.VITE_API_BASE_URL ?? "http://localhost:8000/api").replace(/\/$/, "")
const BASES: Base[] = ["A", "C", "G", "T"]

const IUPAC_CODES: Record<string, Base[]> = {
  A: ["A"],
  C: ["C"],
  G: ["G"],
  T: ["T"],
  R: ["A", "G"],
  Y: ["C", "T"],
  S: ["C", "G"],
  W: ["A", "T"],
  K: ["G", "T"],
  M: ["A", "C"],
  B: ["C", "G", "T"],
  D: ["A", "G", "T"],
  H: ["A", "C", "T"],
  V: ["A", "C", "G"],
  N: ["A", "C", "G", "T"],
}

const BASE_SWATCHES: Record<Base, string> = {
  A: "border-emerald-200 bg-emerald-100 text-emerald-800",
  C: "border-sky-200 bg-sky-100 text-sky-800",
  G: "border-amber-200 bg-amber-100 text-amber-800",
  T: "border-rose-200 bg-rose-100 text-rose-800",
}

const CONTROL_CLASS =
  "h-10 w-full rounded-[8px] border border-stone-300 bg-white px-3 text-sm text-stone-950 shadow-sm outline-none transition focus:border-emerald-600 focus:ring-3 focus:ring-emerald-600/15"

function App() {
  const [file, setFile] = useState<File | null>(null)
  const [profile, setProfile] = useState<FastaProfile | null>(null)
  const [motifPattern, setMotifPattern] = useState("TATAWA")
  const [runMode, setRunMode] = useState<RunMode>("pvalue")
  const [engineMode, setEngineMode] = useState<EngineMode>("custom")
  const [pValueThreshold, setPValueThreshold] = useState("1e-4")
  const [fimoPValueThreshold, setFimoPValueThreshold] = useState("1e-4")
  const [scoreThreshold, setScoreThreshold] = useState("")
  const [scanReverseComplement, setScanReverseComplement] = useState(true)
  const [showAdvanced, setShowAdvanced] = useState(false)
  const [isLoading, setIsLoading] = useState(false)
  const [error, setError] = useState("")
  const [result, setResult] = useState<ScanResult | null>(null)
  const [activeTab, setActiveTab] = useState<ViewTab>("overview")
  const [activeDataset, setActiveDataset] = useState<"custom" | "fimo">("custom")
  const [searchTerm, setSearchTerm] = useState("")
  const [sortMode, setSortMode] = useState<SortMode>("pvalue")
  const [strandFilter, setStrandFilter] = useState<StrandFilter>("all")

  const motifColumns = useMemo(() => parseMotifPattern(motifPattern), [motifPattern])
  const datasetView = useMemo(() => {
    if (!result) return null
    if (result.engine === "both") {
      if (activeDataset === "fimo" && result.fimo) return result.fimo
      if (result.custom) return result.custom
    }
    return { metadata: result.metadata, matches: result.matches }
  }, [activeDataset, result])

  const matches = useMemo(() => datasetView?.matches ?? [], [datasetView])
  const stats = useMemo(() => calculateStats(matches), [matches])
  const geneSummary = useMemo(() => summarizeGenes(matches), [matches])

  const filteredMatches = useMemo(() => {
    const query = searchTerm.trim().toLowerCase()
    const filtered = matches.filter((match) => {
      const passesSearch =
        !query ||
        match.gene.toLowerCase().includes(query) ||
        match.matched_sequence.toLowerCase().includes(query)
      const passesStrand = strandFilter === "all" || match.source_strand === strandFilter
      return passesSearch && passesStrand
    })
    return [...filtered].sort((a, b) => sortMatches(a, b, sortMode))
  }, [matches, searchTerm, sortMode, strandFilter])

  const visibleMatches = filteredMatches.slice(0, 300)
  const currentStep = result ? 3 : isLoading ? 2 : 1

  async function handleFileChange(event: ChangeEvent<HTMLInputElement>) {
    const selectedFile = event.target.files?.[0] ?? null
    setFile(selectedFile)
    setError("")
    setResult(null)

    if (!selectedFile) {
      setProfile(null)
      return
    }

    try {
      const text = await selectedFile.text()
      setProfile(analyzeFasta(text, selectedFile))
    } catch {
      setProfile(null)
      setError("Could not read FASTA file in browser.")
    }
  }

  async function handleSubmit(event: FormEvent<HTMLFormElement>) {
    event.preventDefault()
    setError("")

    if (!file) {
      setError("Attach one FASTA file before running.")
      return
    }

    const trimmedMotif = motifPattern.trim()
    if (!trimmedMotif) {
      setError("Enter a motif pattern.")
      return
    }

    if (!motifColumns.length) {
      setError("Motif pattern is invalid. Use IUPAC symbols and optional groups like [AT].")
      return
    }

    if (runMode === "pvalue") {
      const pValue = Number(pValueThreshold)
      if (!Number.isFinite(pValue) || pValue <= 0 || pValue > 1) {
        setError("P-value threshold must be > 0 and <= 1.")
        return
      }
    } else {
      const score = Number(scoreThreshold)
      if (!Number.isInteger(score)) {
        setError("Score threshold must be an integer.")
        return
      }
    }

    if (engineMode === "fimo" || engineMode === "both") {
      const fimoThreshold = Number(fimoPValueThreshold)
      if (!Number.isFinite(fimoThreshold) || fimoThreshold <= 0 || fimoThreshold > 1) {
        setError("FIMO p-value threshold must be > 0 and <= 1.")
        return
      }
    }

    const formData = new FormData()
    formData.append("file", file)
    formData.append("motif_pattern", trimmedMotif)
    formData.append("scan_reverse_complement", String(scanReverseComplement))
    formData.append("engine", engineMode)
    formData.append("fimo_p_value_threshold", fimoPValueThreshold)
    if (runMode === "pvalue") {
      formData.append("p_value_threshold", pValueThreshold)
    } else {
      formData.append("score_threshold", scoreThreshold)
    }

    setIsLoading(true)
    setResult(null)

    try {
      const response = await fetch(`${API_BASE}/motifs/`, { method: "POST", body: formData })
      const raw = await response.text()
      if (!response.ok) {
        throw new Error(readApiError(raw))
      }

      const payload = parseScanResponse(raw)
      const completedAt = new Date().toISOString()

      if (payload.engine === "both") {
        const customMatches = payload.custom.matches.map((match) => ({
          ...match,
          negLog10PValue: match["-log10_pval"],
        }))
        const fimoMatches = payload.fimo.matches.map((match) => ({
          ...match,
          negLog10PValue: match["-log10_pval"],
        }))

        setResult({
          engine: "both",
          metadata: payload.custom.metadata,
          matches: customMatches,
          custom: {
            metadata: payload.custom.metadata,
            matches: customMatches,
          },
          fimo: {
            metadata: payload.fimo.metadata,
            matches: fimoMatches,
          },
          comparison: payload.comparison,
          completedAt,
        })
        setActiveDataset("custom")
        setActiveTab("compare")
      } else {
        const mappedMatches = payload.matches.map((match) => ({
          ...match,
          negLog10PValue: match["-log10_pval"],
        }))

        setResult({
          engine: payload.engine ?? engineMode,
          metadata: payload.metadata,
          matches: mappedMatches,
          completedAt,
        })
        setActiveDataset(payload.engine === "fimo" ? "fimo" : "custom")
        setActiveTab("overview")
      }
    } catch (caughtError) {
      setError(caughtError instanceof Error ? caughtError.message : "Scan failed.")
    } finally {
      setIsLoading(false)
    }
  }

  function exportJson() {
    if (!result) return
    downloadText(
      `motif-scan-${safeTimestamp()}.json`,
      JSON.stringify({ profile, ...result }, null, 2),
      "application/json",
    )
  }

  function exportCsv() {
    if (!result) return
    const header = [
      "gene",
      "matched_sequence",
      "source_position",
      "source_end_position",
      "source_strand",
      "score",
      "p_value",
      "neg_log10_p",
    ]
    const rows = filteredMatches.map((match) => [
      match.gene,
      match.matched_sequence,
      String(match.source_position),
      String(match.source_end_position),
      match.source_strand,
      String(match.score),
      match.p_value,
      String(match.negLog10PValue ?? ""),
    ])
    downloadText(
      `motif-matches-${safeTimestamp()}.csv`,
      [header, ...rows].map((row) => row.map(escapeCsv).join(",")).join("\n"),
      "text/csv",
    )
  }

  return (
    <main className="min-h-screen bg-[radial-gradient(circle_at_top,_#eefaf3_0%,_#f8fbf7_48%,_#ffffff_100%)] text-stone-950">
      <header className="border-b border-stone-200/80 bg-white/80 backdrop-blur">
        <div className="mx-auto flex max-w-6xl items-center justify-between gap-3 px-4 py-4 sm:px-6">
          <div className="flex items-center gap-3">
            <div className="flex size-10 items-center justify-center rounded-[8px] border border-emerald-200 bg-emerald-50 text-emerald-700">
              <Microscope className="size-5" />
            </div>
            <div>
              <p className="text-xs font-semibold uppercase text-emerald-700">Step by step motif scan</p>
              <h1 className="text-lg font-semibold sm:text-xl">PWM Motif Explorer</h1>
            </div>
          </div>
          <div className="hidden items-center gap-2 rounded-[8px] border border-stone-200 bg-stone-50 px-3 py-2 text-xs text-stone-600 sm:flex">
            <Activity className="size-4 text-emerald-700" />
            <code>{API_BASE}/motifs/</code>
          </div>
        </div>
      </header>

      <div className="mx-auto max-w-6xl px-4 py-6 sm:px-6">
        <StepIndicator currentStep={currentStep} />

        <form className="mt-6 space-y-5" onSubmit={handleSubmit}>
          <section className="rounded-[10px] border border-stone-200 bg-white p-5 shadow-sm">
            <div className="mx-auto max-w-3xl">
              <p className="text-xs font-semibold uppercase text-emerald-700">Step 1</p>
              <h2 className="mt-1 text-xl font-semibold">Enter motif and attach FASTA</h2>
              <p className="mt-1 text-sm text-stone-600">
                Start with one motif input in the center. The base viewer updates while you type.
              </p>

              <div className="mt-5">
                <label className="block text-sm font-medium text-stone-800" htmlFor="motif-pattern">
                  Motif pattern (IUPAC + optional [group])
                </label>
                <input
                  id="motif-pattern"
                  className={`${CONTROL_CLASS} mt-2 text-center font-mono text-base uppercase tracking-[0.08em]`}
                  value={motifPattern}
                  onChange={(event) => setMotifPattern(event.target.value.toUpperCase())}
                  placeholder="TATAWA"
                />
                <MotifPreview columns={motifColumns} motifPattern={motifPattern} />
              </div>

              <div className="mt-5">
                <input
                  className="sr-only"
                  id="fasta-file"
                  type="file"
                  accept=".fasta,.fa,.fna"
                  onChange={handleFileChange}
                />
                <label
                  className="flex min-h-28 cursor-pointer flex-col items-center justify-center gap-2 rounded-[10px] border border-dashed border-emerald-300 bg-emerald-50/60 px-4 py-5 text-center transition hover:border-emerald-500 hover:bg-emerald-50"
                  htmlFor="fasta-file"
                >
                  <UploadCloud className="size-6 text-emerald-700" />
                  <span className="max-w-full truncate text-sm font-semibold text-stone-900">
                    {file ? file.name : "Attach FASTA file"}
                  </span>
                  <span className="text-xs text-stone-600">Accepts .fasta .fa .fna</span>
                </label>
              </div>

              {profile ? <FastaProfilePanel profile={profile} /> : null}

              <div className="mt-5 rounded-[8px] border border-stone-200 bg-stone-50 p-3">
                <button
                  className="flex w-full items-center justify-between text-sm font-medium text-stone-800"
                  type="button"
                  onClick={() => setShowAdvanced((prev) => !prev)}
                >
                  Advanced options
                  <ChevronRight className={`size-4 transition ${showAdvanced ? "rotate-90" : ""}`} />
                </button>

                {showAdvanced ? (
                  <div className="mt-3 space-y-3 border-t border-stone-200 pt-3">
                    <div>
                      <span className="block text-sm font-medium text-stone-800">Scanner engine</span>
                      <div className="mt-2 grid grid-cols-3 gap-1 rounded-[8px] border border-stone-200 bg-stone-100 p-1">
                        <button
                          className={segmentClass(engineMode === "custom")}
                          type="button"
                          onClick={() => setEngineMode("custom")}
                        >
                          Custom
                        </button>
                        <button
                          className={segmentClass(engineMode === "fimo")}
                          type="button"
                          onClick={() => setEngineMode("fimo")}
                        >
                          FIMO
                        </button>
                        <button
                          className={segmentClass(engineMode === "both")}
                          type="button"
                          onClick={() => setEngineMode("both")}
                        >
                          Both
                        </button>
                      </div>
                    </div>

                    <div>
                      <span className="block text-sm font-medium text-stone-800">Threshold mode</span>
                      <div className="mt-2 grid grid-cols-2 gap-1 rounded-[8px] border border-stone-200 bg-stone-100 p-1">
                        <button
                          className={segmentClass(runMode === "pvalue")}
                          type="button"
                          onClick={() => setRunMode("pvalue")}
                        >
                          P-value
                        </button>
                        <button
                          className={segmentClass(runMode === "score")}
                          type="button"
                          onClick={() => setRunMode("score")}
                        >
                          Score
                        </button>
                      </div>
                    </div>

                    {runMode === "pvalue" ? (
                      <div>
                        <label className="block text-sm font-medium text-stone-800" htmlFor="pvalue-input">
                          P-value threshold
                        </label>
                        <input
                          id="pvalue-input"
                          className={`${CONTROL_CLASS} mt-2 font-mono`}
                          value={pValueThreshold}
                          onChange={(event) => setPValueThreshold(event.target.value)}
                        />
                      </div>
                    ) : (
                      <div>
                        <label className="block text-sm font-medium text-stone-800" htmlFor="score-input">
                          Integer score threshold
                        </label>
                        <input
                          id="score-input"
                          className={`${CONTROL_CLASS} mt-2 font-mono`}
                          value={scoreThreshold}
                          onChange={(event) => setScoreThreshold(event.target.value)}
                          placeholder="120"
                        />
                      </div>
                    )}

                    <label className="flex items-center justify-between rounded-[8px] border border-stone-200 bg-white px-3 py-2">
                      <span className="text-sm font-medium text-stone-800">Scan reverse complement</span>
                      <input
                        className="size-5 accent-emerald-700"
                        type="checkbox"
                        checked={scanReverseComplement}
                        onChange={(event) => setScanReverseComplement(event.target.checked)}
                      />
                    </label>

                    {(engineMode === "fimo" || engineMode === "both") ? (
                      <div>
                        <label className="block text-sm font-medium text-stone-800" htmlFor="fimo-pvalue-input">
                          FIMO p-value threshold
                        </label>
                        <input
                          id="fimo-pvalue-input"
                          className={`${CONTROL_CLASS} mt-2 font-mono`}
                          value={fimoPValueThreshold}
                          onChange={(event) => setFimoPValueThreshold(event.target.value)}
                        />
                      </div>
                    ) : null}
                  </div>
                ) : null}
              </div>

              <Button className="mt-5 h-10 w-full bg-emerald-700 text-white hover:bg-emerald-800" disabled={isLoading}>
                {isLoading ? <Loader2 className="size-4 animate-spin" /> : <Dna className="size-4" />}
                {isLoading ? "Processing..." : "Start processing"}
              </Button>
            </div>
          </section>
        </form>

        {isLoading ? (
          <section className="mt-5 rounded-[10px] border border-stone-200 bg-white p-5 shadow-sm">
            <p className="text-xs font-semibold uppercase text-emerald-700">Step 2</p>
            <div className="mt-2 flex items-center gap-3 text-sm text-stone-700">
              <Loader2 className="size-4 animate-spin text-emerald-700" />
              Running PWM scan with current motif and FASTA.
            </div>
          </section>
        ) : null}

        {error ? (
          <section className="mt-5 rounded-[10px] border border-rose-200 bg-rose-50 p-4 text-sm text-rose-900">
            {error}
          </section>
        ) : null}

        {result ? (
          <section className="mt-5 rounded-[10px] border border-stone-200 bg-white shadow-sm">
            <div className="border-b border-stone-200 p-4">
              <p className="text-xs font-semibold uppercase text-emerald-700">Step 3</p>
              <div className="mt-1 flex flex-wrap items-center justify-between gap-3">
                <h2 className="text-xl font-semibold">Results</h2>
                <div className="flex gap-2">
                  <Button className="bg-white" type="button" variant="outline" onClick={exportCsv}>
                    <FileSpreadsheet className="size-4" />
                    CSV
                  </Button>
                  <Button className="bg-white" type="button" variant="outline" onClick={exportJson}>
                    <FileJson className="size-4" />
                    JSON
                  </Button>
                </div>
              </div>
            </div>

            <div className="p-4">
              <div className="grid gap-3 sm:grid-cols-2 lg:grid-cols-4">
                <MetricCard label="Total matches" value={formatInteger(stats.totalMatches)} />
                <MetricCard label="Genes with hits" value={formatInteger(stats.uniqueGenes)} />
                <MetricCard label="Best p-value" value={stats.bestPValue === null ? "-" : formatPValue(stats.bestPValue)} />
                <MetricCard label="Max score" value={stats.maxScore === null ? "-" : String(stats.maxScore)} />
              </div>

              <div className="mt-4 flex flex-wrap items-center gap-2">
                <TabButton active={activeTab === "overview"} onClick={() => setActiveTab("overview")}>
                  Overview
                </TabButton>
                <TabButton active={activeTab === "genes"} onClick={() => setActiveTab("genes")}>
                  Genes
                </TabButton>
                <TabButton active={activeTab === "matches"} onClick={() => setActiveTab("matches")}>
                  Matches
                </TabButton>
                {result.engine === "both" ? (
                  <TabButton active={activeTab === "compare"} onClick={() => setActiveTab("compare")}>
                    Compare
                  </TabButton>
                ) : null}
              </div>

              {result.engine === "both" ? (
                <div className="mt-3 inline-flex rounded-[8px] border border-stone-200 bg-stone-100 p-1">
                  <button
                    className={segmentClass(activeDataset === "custom")}
                    onClick={() => setActiveDataset("custom")}
                    type="button"
                  >
                    Custom dataset
                  </button>
                  <button
                    className={segmentClass(activeDataset === "fimo")}
                    onClick={() => setActiveDataset("fimo")}
                    type="button"
                  >
                    FIMO dataset
                  </button>
                </div>
              ) : null}

              {activeTab === "overview" && datasetView ? (
                <OverviewPanel
                  metadata={datasetView.metadata}
                  profile={profile}
                  stats={stats}
                  completedAt={result.completedAt}
                />
              ) : null}

              {activeTab === "genes" ? <GenesPanel rows={geneSummary} /> : null}

              {activeTab === "matches" ? (
                <MatchesPanel
                  searchTerm={searchTerm}
                  setSearchTerm={setSearchTerm}
                  sortMode={sortMode}
                  setSortMode={setSortMode}
                  strandFilter={strandFilter}
                  setStrandFilter={setStrandFilter}
                  filteredMatches={filteredMatches}
                  visibleMatches={visibleMatches}
                />
              ) : null}

              {activeTab === "compare" && result.engine === "both" && result.custom && result.fimo && result.comparison ? (
                <ComparePanel
                  comparison={result.comparison}
                  customMatches={result.custom.matches}
                  fimoMatches={result.fimo.matches}
                />
              ) : null}
            </div>
          </section>
        ) : null}
      </div>
    </main>
  )
}

function StepIndicator({ currentStep }: { currentStep: 1 | 2 | 3 }) {
  return (
    <div className="grid gap-2 rounded-[10px] border border-stone-200 bg-white p-3 shadow-sm sm:grid-cols-3">
      <StepNode index={1} active={currentStep === 1} done={currentStep > 1} label="Setup" />
      <StepNode index={2} active={currentStep === 2} done={currentStep > 2} label="Processing" />
      <StepNode index={3} active={currentStep === 3} done={false} label="Results" />
    </div>
  )
}

function StepNode({ index, active, done, label }: { index: number; active: boolean; done: boolean; label: string }) {
  return (
    <div
      className={`flex items-center gap-2 rounded-[8px] border px-3 py-2 ${
        active
          ? "border-emerald-300 bg-emerald-50"
          : done
            ? "border-emerald-200 bg-white"
            : "border-stone-200 bg-white"
      }`}
    >
      <span
        className={`flex size-6 items-center justify-center rounded-full text-xs font-semibold ${
          done ? "bg-emerald-600 text-white" : active ? "bg-emerald-700 text-white" : "bg-stone-200 text-stone-700"
        }`}
      >
        {done ? <CheckCircle2 className="size-3.5" /> : index}
      </span>
      <span className={`text-sm font-medium ${active ? "text-emerald-800" : "text-stone-700"}`}>{label}</span>
    </div>
  )
}

function MetricCard({ label, value }: { label: string; value: string }) {
  return (
    <div className="rounded-[8px] border border-stone-200 bg-stone-50 p-3">
      <p className="text-xs font-medium uppercase text-stone-500">{label}</p>
      <p className="mt-1 text-xl font-semibold text-stone-900">{value}</p>
    </div>
  )
}

function FastaProfilePanel({ profile }: { profile: FastaProfile }) {
  return (
    <div className="mt-4 rounded-[8px] border border-stone-200 bg-stone-50 p-3">
      <p className="text-sm font-semibold text-stone-900">FASTA snapshot</p>
      <dl className="mt-2 grid gap-2 text-sm sm:grid-cols-2">
        <StatLine label="File size" value={formatBytes(profile.sizeBytes)} />
        <StatLine label="Records" value={formatInteger(profile.records)} />
        <StatLine label="Total bases" value={formatInteger(profile.totalBases)} />
        <StatLine label="Ambiguous bases" value={formatInteger(profile.ambiguousBases)} />
        <StatLine label="GC ratio" value={formatPercent(profile.gcRatio)} />
        <StatLine label="First header" value={profile.firstHeader || "-"} />
      </dl>
    </div>
  )
}

function MotifPreview({ columns, motifPattern }: { columns: Base[][]; motifPattern: string }) {
  if (!motifPattern.trim()) {
    return <p className="mt-3 text-sm text-stone-500">Type a motif to see base composition.</p>
  }

  if (!columns.length) {
    return <p className="mt-3 text-sm text-rose-700">Invalid motif pattern.</p>
  }

  return (
    <div className="mt-3 overflow-x-auto rounded-[8px] border border-stone-200 bg-stone-50 p-3">
      <div className="flex min-w-fit gap-2">
        {columns.map((column, index) => (
          <div className="w-[66px] rounded-[8px] border border-stone-200 bg-white p-2" key={`${index}-${column.join("")}`}>
            <p className="text-center text-[10px] font-semibold uppercase text-stone-500">pos {index + 1}</p>
            <div className="mt-2 space-y-1">
              {BASES.map((base) => {
                const probability = column.includes(base) ? 1 / column.length : 0
                return (
                  <div className="grid grid-cols-[14px_minmax(0,1fr)] items-center gap-1" key={base}>
                    <span className={`rounded-[4px] border text-center text-[10px] font-bold ${BASE_SWATCHES[base]}`}>
                      {base}
                    </span>
                    <div className="h-2 rounded-full bg-stone-100">
                      <div className="h-2 rounded-full bg-emerald-600" style={{ width: `${probability * 100}%` }} />
                    </div>
                  </div>
                )
              })}
            </div>
          </div>
        ))}
      </div>
    </div>
  )
}

function OverviewPanel({
  metadata,
  profile,
  stats,
  completedAt,
}: {
  metadata: MotifMetadata
  profile: FastaProfile | null
  stats: RunStats
  completedAt: string
}) {
  return (
    <div className="mt-4 grid gap-4 lg:grid-cols-[minmax(0,1fr)_340px]">
      <div className="rounded-[8px] border border-stone-200 bg-white p-4">
        <h3 className="text-sm font-semibold text-stone-900">Run summary</h3>
        <dl className="mt-3 grid gap-2 text-sm sm:grid-cols-2">
          <StatLine label="Input file" value={metadata.filename} />
          <StatLine label="Motif pattern" value={metadata.motif_pattern} />
          <StatLine label="Motif length" value={String(metadata.motif_length)} />
          <StatLine label="Cutoff mode" value={metadata.cutoff_mode} />
          <StatLine
            label="Score threshold"
            value={metadata.score_threshold === null ? "derived from p-value" : String(metadata.score_threshold)}
          />
          <StatLine label="Completed at" value={new Date(completedAt).toLocaleString()} />
          <StatLine label="Forward strand hits" value={formatInteger(stats.forwardSource)} />
          <StatLine label="Reverse strand hits" value={formatInteger(stats.reverseSource)} />
          <StatLine label="Canonical bases in file" value={profile ? formatInteger(profile.canonicalBases) : "-"} />
        </dl>
      </div>

      <div className="rounded-[8px] border border-stone-200 bg-stone-50 p-4">
        <h3 className="text-sm font-semibold text-stone-900">Background frequencies</h3>
        <div className="mt-3 space-y-3">
          {BASES.map((base) => (
            <div className="grid grid-cols-[24px_minmax(0,1fr)_50px] items-center gap-2" key={base}>
              <span className={`flex size-6 items-center justify-center rounded-[6px] border text-xs font-bold ${BASE_SWATCHES[base]}`}>
                {base}
              </span>
              <div className="h-2 rounded-full bg-stone-200">
                <div
                  className="h-2 rounded-full bg-emerald-700"
                  style={{ width: `${Math.max(2, (metadata.background_frequencies[base] ?? 0) * 100)}%` }}
                />
              </div>
              <span className="text-right font-mono text-xs text-stone-700">
                {formatPercent(metadata.background_frequencies[base] ?? 0)}
              </span>
            </div>
          ))}
        </div>
      </div>
    </div>
  )
}

function GenesPanel({ rows }: { rows: GeneSummaryRow[] }) {
  return (
    <div className="mt-4 overflow-x-auto rounded-[8px] border border-stone-200">
      <table className="min-w-[720px] w-full text-left text-sm">
        <thead className="bg-stone-100 text-xs uppercase text-stone-600">
          <tr>
            <th className="px-3 py-3">Gene</th>
            <th className="px-3 py-3">Hits</th>
            <th className="px-3 py-3">Best p-value</th>
            <th className="px-3 py-3">Best score</th>
            <th className="px-3 py-3">First position</th>
            <th className="px-3 py-3">Strands</th>
          </tr>
        </thead>
        <tbody className="divide-y divide-stone-200 bg-white">
          {rows.slice(0, 300).map((row) => (
            <tr className="hover:bg-emerald-50/50" key={row.gene}>
              <td className="max-w-80 truncate px-3 py-3 font-medium text-stone-900">{row.gene}</td>
              <td className="px-3 py-3 font-mono">{row.hits}</td>
              <td className="px-3 py-3 font-mono">{formatPValue(row.bestPValue)}</td>
              <td className="px-3 py-3 font-mono">{row.bestScore}</td>
              <td className="px-3 py-3 font-mono">{row.firstPosition}</td>
              <td className="px-3 py-3 text-stone-700">{row.sourceStrands}</td>
            </tr>
          ))}
          {!rows.length ? (
            <tr>
              <td className="px-3 py-12 text-center text-sm text-stone-500" colSpan={6}>
                No genes to summarize.
              </td>
            </tr>
          ) : null}
        </tbody>
      </table>
    </div>
  )
}

function ComparePanel({
  comparison,
  customMatches,
  fimoMatches,
}: {
  comparison: ComparisonStats
  customMatches: MotifMatch[]
  fimoMatches: MotifMatch[]
}) {
  const customBest = customMatches.length
    ? [...customMatches].sort((a, b) => a.p_value_numeric - b.p_value_numeric)[0]
    : null
  const fimoBest = fimoMatches.length
    ? [...fimoMatches].sort((a, b) => a.p_value_numeric - b.p_value_numeric)[0]
    : null

  return (
    <div className="mt-4 grid gap-4 lg:grid-cols-2">
      <div className="rounded-[8px] border border-stone-200 bg-white p-4">
        <h3 className="text-sm font-semibold text-stone-900">Overlap summary</h3>
        <dl className="mt-3 space-y-2 text-sm">
          <StatLine label="Custom total" value={formatInteger(comparison.custom_total)} />
          <StatLine label="FIMO total" value={formatInteger(comparison.fimo_total)} />
          <StatLine label="Shared hits" value={formatInteger(comparison.overlap_total)} />
          <StatLine label="Only custom" value={formatInteger(comparison.only_custom_total)} />
          <StatLine label="Only FIMO" value={formatInteger(comparison.only_fimo_total)} />
          <StatLine label="Custom overlap ratio" value={formatPercent(comparison.overlap_ratio_custom)} />
          <StatLine label="FIMO overlap ratio" value={formatPercent(comparison.overlap_ratio_fimo)} />
        </dl>
      </div>

      <div className="rounded-[8px] border border-stone-200 bg-stone-50 p-4">
        <h3 className="text-sm font-semibold text-stone-900">Best hit by engine</h3>
        <div className="mt-3 grid gap-3">
          <EngineBestHit
            title="Custom"
            gene={customBest?.gene ?? "-"}
            pValue={customBest ? formatPValue(customBest.p_value_numeric) : "-"}
            score={customBest ? String(customBest.score) : "-"}
          />
          <EngineBestHit
            title="FIMO"
            gene={fimoBest?.gene ?? "-"}
            pValue={fimoBest ? formatPValue(fimoBest.p_value_numeric) : "-"}
            score={fimoBest ? String(fimoBest.score) : "-"}
          />
        </div>
      </div>
    </div>
  )
}

function EngineBestHit({
  title,
  gene,
  pValue,
  score,
}: {
  title: string
  gene: string
  pValue: string
  score: string
}) {
  return (
    <div className="rounded-[8px] border border-stone-200 bg-white p-3">
      <p className="text-xs font-semibold uppercase text-stone-600">{title}</p>
      <dl className="mt-2 space-y-1 text-sm">
        <StatLine label="Gene" value={gene} />
        <StatLine label="Best p-value" value={pValue} />
        <StatLine label="Score" value={score} />
      </dl>
    </div>
  )
}

function MatchesPanel({
  searchTerm,
  setSearchTerm,
  sortMode,
  setSortMode,
  strandFilter,
  setStrandFilter,
  filteredMatches,
  visibleMatches,
}: {
  searchTerm: string
  setSearchTerm: (value: string) => void
  sortMode: SortMode
  setSortMode: (value: SortMode) => void
  strandFilter: StrandFilter
  setStrandFilter: (value: StrandFilter) => void
  filteredMatches: MotifMatch[]
  visibleMatches: MotifMatch[]
}) {
  return (
    <div className="mt-4">
      <div className="grid gap-3 md:grid-cols-[minmax(0,1fr)_170px_170px]">
        <label className="relative block">
          <Search className="pointer-events-none absolute left-3 top-1/2 size-4 -translate-y-1/2 text-stone-500" />
          <input
            className={`${CONTROL_CLASS} pl-9`}
            placeholder="Search gene or sequence"
            value={searchTerm}
            onChange={(event) => setSearchTerm(event.target.value)}
          />
        </label>
        <label className="relative block">
          <Filter className="pointer-events-none absolute left-3 top-1/2 size-4 -translate-y-1/2 text-stone-500" />
          <select
            className={`${CONTROL_CLASS} appearance-none pl-9`}
            value={strandFilter}
            onChange={(event) => setStrandFilter(event.target.value as StrandFilter)}
          >
            <option value="all">All strands</option>
            <option value="forward">Forward source</option>
            <option value="reverse">Reverse source</option>
          </select>
        </label>
        <label className="relative block">
          <ArrowDownUp className="pointer-events-none absolute left-3 top-1/2 size-4 -translate-y-1/2 text-stone-500" />
          <select
            className={`${CONTROL_CLASS} appearance-none pl-9`}
            value={sortMode}
            onChange={(event) => setSortMode(event.target.value as SortMode)}
          >
            <option value="pvalue">P-value</option>
            <option value="score">Score</option>
            <option value="gene">Gene</option>
            <option value="position">Position</option>
          </select>
        </label>
      </div>

      <div className="mt-3 overflow-x-auto rounded-[8px] border border-stone-200">
        <table className="min-w-[940px] w-full text-left text-sm">
          <thead className="bg-stone-100 text-xs uppercase text-stone-600">
            <tr>
              <th className="px-3 py-3">Gene</th>
              <th className="px-3 py-3">Sequence</th>
              <th className="px-3 py-3">Coordinates</th>
              <th className="px-3 py-3">Strand</th>
              <th className="px-3 py-3">Score</th>
              <th className="px-3 py-3">P-value</th>
              <th className="px-3 py-3">-log10 p</th>
            </tr>
          </thead>
          <tbody className="divide-y divide-stone-200 bg-white">
            {visibleMatches.map((match, index) => (
              <tr className="hover:bg-emerald-50/50" key={`${match.gene}-${match.source_position}-${index}`}>
                <td className="max-w-56 truncate px-3 py-3 font-medium text-stone-900">{match.gene}</td>
                <td className="px-3 py-3">
                  <SequenceSwatches sequence={match.matched_sequence} />
                </td>
                <td className="px-3 py-3 font-mono text-xs text-stone-700">
                  {match.source_position}-{match.source_end_position}
                </td>
                <td className="px-3 py-3">
                  <span
                    className={`inline-flex items-center rounded-[8px] border px-2 py-1 text-xs font-medium ${
                      match.source_strand === "reverse"
                        ? "border-rose-200 bg-rose-50 text-rose-800"
                        : "border-emerald-200 bg-emerald-50 text-emerald-800"
                    }`}
                  >
                    {match.source_strand_symbol} {match.source_strand}
                  </span>
                </td>
                <td className="px-3 py-3 font-mono">{match.score}</td>
                <td className="px-3 py-3 font-mono">{match.p_value}</td>
                <td className="px-3 py-3 font-mono">{match.negLog10PValue?.toFixed(2) ?? "-"}</td>
              </tr>
            ))}
            {!visibleMatches.length ? (
              <tr>
                <td className="px-3 py-12 text-center text-sm text-stone-500" colSpan={7}>
                  No motif hits for current filters.
                </td>
              </tr>
            ) : null}
          </tbody>
        </table>
      </div>

      <p className="mt-2 text-xs text-stone-500">
        Showing {formatInteger(visibleMatches.length)} of {formatInteger(filteredMatches.length)} matches.
      </p>
    </div>
  )
}

function TabButton({
  active,
  children,
  onClick,
}: {
  active: boolean
  children: string
  onClick: () => void
}) {
  return (
    <button
      className={`h-8 rounded-[8px] border px-3 text-sm font-medium transition ${
        active
          ? "border-emerald-300 bg-emerald-50 text-emerald-800"
          : "border-stone-200 bg-white text-stone-700 hover:bg-stone-50"
      }`}
      onClick={onClick}
      type="button"
    >
      {children}
    </button>
  )
}

function StatLine({ label, value }: { label: string; value: string }) {
  return (
    <div className="flex items-center justify-between gap-3">
      <dt className="truncate text-stone-500">{label}</dt>
      <dd className="truncate text-right font-mono text-stone-900">{value}</dd>
    </div>
  )
}

function SequenceSwatches({ sequence }: { sequence: string }) {
  return (
    <div className="flex flex-wrap gap-1 font-mono">
      {sequence.split("").map((base, index) => {
        const typedBase = BASES.includes(base as Base) ? (base as Base) : null
        return (
          <span
            className={`flex size-6 items-center justify-center rounded-[6px] border text-xs font-bold ${
              typedBase ? BASE_SWATCHES[typedBase] : "border-stone-200 bg-stone-100 text-stone-600"
            }`}
            key={`${base}-${index}`}
          >
            {base}
          </span>
        )
      })}
    </div>
  )
}

function parseMotifPattern(pattern: string): Base[][] {
  const columns: Base[][] = []
  const normalized = pattern.trim().toUpperCase()
  let index = 0

  while (index < normalized.length) {
    const char = normalized[index]

    if (char === "[") {
      const end = normalized.indexOf("]", index)
      if (end === -1) return []

      const group = normalized.slice(index + 1, end)
      const bases = [...new Set(group.split(""))].filter((base): base is Base => BASES.includes(base as Base))
      if (!bases.length || bases.length !== group.length) return []

      columns.push(bases)
      index = end + 1
      continue
    }

    const bases = IUPAC_CODES[char]
    if (!bases) return []
    columns.push(bases)
    index += 1
  }

  return columns
}

function analyzeFasta(text: string, file: File): FastaProfile {
  const counts: Record<Base, number> = { A: 0, C: 0, G: 0, T: 0 }
  let records = 0
  let ambiguousBases = 0
  let firstHeader = ""
  let sawSequence = false

  for (const rawLine of text.split(/\r?\n/)) {
    const line = rawLine.trim()
    if (!line) continue

    if (line.startsWith(">")) {
      records += 1
      if (!firstHeader) firstHeader = line
      continue
    }

    sawSequence = true
    for (const char of line.toUpperCase()) {
      if (char === "A" || char === "C" || char === "G" || char === "T") {
        counts[char] += 1
      } else if (/[A-Z]/.test(char)) {
        ambiguousBases += 1
      }
    }
  }

  const canonicalBases = BASES.reduce((sum, base) => sum + counts[base], 0)
  const totalBases = canonicalBases + ambiguousBases
  const inferredRecords = records || (sawSequence ? 1 : 0)
  const gcRatio = canonicalBases ? (counts.G + counts.C) / canonicalBases : 0

  return {
    fileName: file.name,
    sizeBytes: file.size,
    records: inferredRecords,
    totalBases,
    canonicalBases,
    ambiguousBases,
    gcRatio,
    counts,
    firstHeader,
  }
}

function calculateStats(matches: MotifMatch[]): RunStats {
  const pValues = matches
    .map((match) => match.p_value_numeric)
    .filter((value) => Number.isFinite(value))
    .sort((a, b) => a - b)
  const scores = matches.map((match) => match.score)

  return {
    totalMatches: matches.length,
    uniqueGenes: new Set(matches.map((match) => match.gene)).size,
    forwardSource: matches.filter((match) => match.source_strand === "forward").length,
    reverseSource: matches.filter((match) => match.source_strand === "reverse").length,
    bestPValue: pValues[0] ?? null,
    maxScore: scores.length ? Math.max(...scores) : null,
  }
}

function summarizeGenes(matches: MotifMatch[]): GeneSummaryRow[] {
  const grouped = new Map<string, MotifMatch[]>()

  for (const match of matches) {
    grouped.set(match.gene, [...(grouped.get(match.gene) ?? []), match])
  }

  return [...grouped.entries()]
    .map(([gene, geneMatches]) => {
      const sortedByP = [...geneMatches].sort((a, b) => a.p_value_numeric - b.p_value_numeric)
      const firstByPosition = [...geneMatches].sort((a, b) => a.source_position - b.source_position)[0]
      const sourceStrands = [...new Set(geneMatches.map((match) => match.source_strand_symbol))].join(" / ")
      return {
        gene,
        hits: geneMatches.length,
        bestPValue: sortedByP[0].p_value_numeric,
        bestScore: Math.max(...geneMatches.map((match) => match.score)),
        firstPosition: firstByPosition.source_position,
        sourceStrands,
      }
    })
    .sort((a, b) => a.bestPValue - b.bestPValue || b.hits - a.hits)
}

function sortMatches(a: MotifMatch, b: MotifMatch, sortMode: SortMode) {
  if (sortMode === "score") return b.score - a.score || a.p_value_numeric - b.p_value_numeric
  if (sortMode === "gene") return a.gene.localeCompare(b.gene) || a.source_position - b.source_position
  if (sortMode === "position") return a.source_position - b.source_position || a.gene.localeCompare(b.gene)
  return a.p_value_numeric - b.p_value_numeric || b.score - a.score
}

function parseScanResponse(raw: string): SingleEnginePayload | BothEnginesPayload {
  const trimmed = raw.trim()
  if (!trimmed) throw new Error("API returned an empty response.")

  try {
    return JSON.parse(trimmed) as SingleEnginePayload | BothEnginesPayload
  } catch {
    const payloads = trimmed
      .split(/\r?\n/)
      .filter(Boolean)
      .map((line) => JSON.parse(line) as SingleEnginePayload | BothEnginesPayload)
    const payload = payloads.find((candidate) =>
      "engine" in candidate
        ? candidate.engine === "both"
          ? Boolean(candidate.custom && candidate.fimo)
          : Boolean((candidate as SingleEnginePayload).metadata && (candidate as SingleEnginePayload).matches)
        : Boolean((candidate as SingleEnginePayload).metadata && (candidate as SingleEnginePayload).matches),
    )
    if (!payload) throw new Error("Could not parse API response.")
    return payload
  }
}

function readApiError(raw: string) {
  try {
    const parsed = JSON.parse(raw) as { detail?: string }
    return parsed.detail ?? raw
  } catch {
    return raw || "API rejected scan request."
  }
}

function segmentClass(active: boolean) {
  return `h-8 rounded-[6px] text-sm font-medium transition ${
    active ? "bg-white text-emerald-800 shadow-sm" : "text-stone-600 hover:text-stone-950"
  }`
}

function formatBytes(bytes: number) {
  if (bytes < 1024) return `${bytes} B`
  if (bytes < 1024 * 1024) return `${(bytes / 1024).toFixed(1)} KB`
  return `${(bytes / (1024 * 1024)).toFixed(1)} MB`
}

function formatInteger(value: number) {
  return new Intl.NumberFormat().format(value)
}

function formatPercent(value: number) {
  return `${(value * 100).toFixed(1)}%`
}

function formatPValue(value: number) {
  if (!Number.isFinite(value)) return "-"
  return value.toExponential(2)
}

function safeTimestamp() {
  return new Date().toISOString().replaceAll(":", "-").replace(/\.\d{3}Z$/, "Z")
}

function escapeCsv(value: string) {
  return `"${value.replaceAll('"', '""')}"`
}

function downloadText(filename: string, text: string, type: string) {
  const blob = new Blob([text], { type })
  const url = URL.createObjectURL(blob)
  const anchor = document.createElement("a")
  anchor.href = url
  anchor.download = filename
  anchor.click()
  URL.revokeObjectURL(url)
}

export default App
