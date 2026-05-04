import {ArrowDownUp, Dna, FileJson, FileSpreadsheet, Filter, Loader2, Search, UploadCloud} from "lucide-react"
import type {ChangeEvent, FormEvent} from "react"
import {useMemo, useState} from "react"
import {toast} from "sonner"

import {Button} from "@/components/ui/button"
import {Input} from "@/components/ui/input.tsx";
import {Toaster} from "@/components/ui/sonner"

type Base = "A" | "C" | "G" | "T"
type ViewTab = "overview" | "genes" | "matches"
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

type PValueEngineDiagnostics = {
  method?: string
  initial_granularity?: number
  decreasing_step?: number
  levels_used?: number
  levels?: Array<{
    epsilon: number
    max_error: number
    score_count: number
  }>
}

type MotifMetadata = {
  filename: string
  motif_pattern: string
  motif_length: number
  score_threshold: number | null
  cutoff_mode: string
  background_frequencies: Record<Base, number>
  scan_reverse_complement?: boolean
  custom_pseudo_fraction?: number
  pvalue_engine?: PValueEngineDiagnostics
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
  score_log_odds?: number
  p_value: string
  p_value_numeric: number
  "-log10_pval": number | null
}

type SingleEnginePayload = {
  engine?: "custom"
  metadata: MotifMetadata
  matches: RawMotifMatch[]
}

type ScanStreamEvent =
  | { type: "progress"; progress: number; stage: string }
  | { type: "result"; payload: SingleEnginePayload }
  | { type: "error"; message: string }

type MotifMatch = RawMotifMatch & {
  negLog10PValue: number | null
}

type ScanResult = {
  metadata: MotifMetadata
  matches: MotifMatch[]
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
const DEFAULT_P_VALUE_THRESHOLD = "1e-4"

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

const MOTIF_INPUT_ALLOWED_CHARS = new Set([
  ...Object.keys(IUPAC_CODES),
  "[",
  "]",
])

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
  const [motifPattern, setMotifPattern] = useState("")
  const [isLoading, setIsLoading] = useState(false)
  const [loadingProgress, setLoadingProgress] = useState(0)
  const [loadingStage, setLoadingStage] = useState("Preparing scan")
  const [result, setResult] = useState<ScanResult | null>(null)
  const [activeTab, setActiveTab] = useState<ViewTab>("overview")
  const [searchTerm, setSearchTerm] = useState("")
  const [sortMode, setSortMode] = useState<SortMode>("pvalue")
  const [strandFilter, setStrandFilter] = useState<StrandFilter>("all")

  const motifColumns = useMemo(() => parseMotifPattern(motifPattern), [motifPattern])
  const motifBaseCount = motifColumns.length
  const matches = useMemo(() => result?.matches ?? [], [result])
  const stats = useMemo(() => calculateStats(matches), [matches])
  const geneSummary = useMemo(() => summarizeGenes(matches), [matches])

  const filteredMatches = useMemo(() => {
    const query = searchTerm.trim().toLowerCase()
    const filtered = matches.filter((match) => {
      const passesSearch =
        !query || match.gene.toLowerCase().includes(query) || match.matched_sequence.toLowerCase().includes(query)
      const passesStrand = strandFilter === "all" || match.source_strand === strandFilter
      return passesSearch && passesStrand
    })
    return [...filtered].sort((a, b) => sortMatches(a, b, sortMode))
  }, [matches, searchTerm, sortMode, strandFilter])

  const visibleMatches = filteredMatches.slice(0, 300)

  async function handleFileChange(event: ChangeEvent<HTMLInputElement>) {
    const selectedFile = event.target.files?.[0] ?? null
    setFile(selectedFile)
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
      toast.error("Could not read FASTA file in browser.")
    }
  }

  async function handleSubmit(event: FormEvent<HTMLFormElement>) {
    event.preventDefault()
    toast.dismiss()

    const trimmedMotif = motifPattern.trim()
    if (!trimmedMotif) {
      toast.warning("Enter a motif pattern.")
      return
    }

    if (!motifColumns.length) {
      toast.error("Motif pattern is invalid. Use IUPAC symbols and optional groups like [AT].")
      return
    }

    if (!file) {
      toast.warning("Attach one FASTA file before running.")
      return
    }

    const pValue = Number(DEFAULT_P_VALUE_THRESHOLD)
    if (!Number.isFinite(pValue) || pValue <= 0 || pValue > 1) {
      toast.error("P-value threshold must be > 0 and <= 1.")
      return
    }

    const formData = new FormData()
    formData.append("file", file)
    formData.append("motif_pattern", trimmedMotif)
    formData.append("scan_reverse_complement", "true")
    formData.append("p_value_threshold", DEFAULT_P_VALUE_THRESHOLD)

    setLoadingProgress(0)
    setLoadingStage("Uploading request")
    setIsLoading(true)
    setResult(null)

    let scanCompleted = false
    try {
      const response = await fetch(`${API_BASE}/motifs/stream`, {method: "POST", body: formData})
      if (!response.ok) {
        const raw = await response.text()
        throw new Error(readApiError(raw))
      }

      const payload = await readScanStream(response, (progress, stage) => {
        setLoadingProgress(progress)
        setLoadingStage(stage)
      })
      const mappedMatches = payload.matches.map((match) => ({
        ...match,
        negLog10PValue: match["-log10_pval"],
      }))

      setResult({
        metadata: payload.metadata,
        matches: mappedMatches,
        completedAt: new Date().toISOString(),
      })
      setActiveTab("overview")
      setSearchTerm("")
      setSortMode("pvalue")
      setStrandFilter("all")
      scanCompleted = true
    } catch (caughtError) {
      setLoadingStage("Scan failed")
      toast.error(caughtError instanceof Error ? caughtError.message : "Scan failed.")
    } finally {
      if (scanCompleted) {
        setLoadingProgress(100)
        setLoadingStage("Scan complete")
        await new Promise((resolve) => window.setTimeout(resolve, 180))
      }
      setIsLoading(false)
    }
  }

  function handleMotifPatternChange(value: string) {
    const normalized = value.toUpperCase()
    const sanitized = Array.from(normalized)
      .filter((char) => MOTIF_INPUT_ALLOWED_CHARS.has(char))
      .join("")
    setMotifPattern(sanitized)
  }

  function exportJson() {
    if (!result) return
    downloadText(
      `motif-scan-${safeTimestamp()}.json`,
      JSON.stringify({profile, ...result}, null, 2),
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
      `motif-matches-custom-${safeTimestamp()}.csv`,
      [header, ...rows].map((row) => row.map(escapeCsv).join(",")).join("\n"),
      "text/csv",
    )
  }

  return (
    <main className="min-h-screen bg-emerald-50/50 text-stone-950">
      <Toaster closeButton position="top-center" richColors/>
      <div className="relative z-10 mx-auto max-w-6xl px-4 py-6 sm:px-6">
        {!result && !isLoading ? (
          <section className="flex min-h-[78vh] items-center justify-center">
            <form
              className="w-full max-w-xl"
              onSubmit={handleSubmit}
            >
              <div className="space-y-3">
                <div className="grid gap-2 sm:grid-cols-[minmax(0,1fr)_132px]">
                  <label className="relative block min-w-0" htmlFor="motif-pattern">
                    <Input
                      id="motif-pattern"
                      value={motifPattern}
                      onChange={(event) => handleMotifPatternChange(event.target.value)}
                      placeholder="Input motif pattern"
                      className="
                        h-12 rounded-[8px] border-emerald-100 bg-white px-4 pr-24
                        font-mono text-base font-medium text-stone-950 shadow-sm shadow-stone-900/5
                        placeholder:font-sans placeholder:font-normal placeholder:text-stone-400
                        focus-visible:border-emerald-600
                        focus-visible:ring-3 focus-visible:ring-emerald-600/15
                      "
                    />
                    <span className="pointer-events-none absolute right-3 top-1/2 -translate-y-1/2 rounded-[6px] bg-emerald-50 px-2 py-1 font-mono text-xs font-medium text-emerald-800">
                      {motifBaseCount} {motifBaseCount === 1 ? "base" : "bases"}
                    </span>
                  </label>
                  <Button
                    className="
                      h-12 rounded-[8px] border-none bg-emerald-700 px-4 text-white shadow-sm
                      shadow-emerald-900/15 hover:bg-emerald-800
                    "
                    type="submit"
                  >
                    <Dna className="size-4"/>
                    Run Scan
                  </Button>
                </div>
                <MotifPreview columns={motifColumns} motifPattern={motifPattern}/>

                <input
                  className="sr-only"
                  id="fasta-file"
                  type="file"
                  accept=".fasta,.fa,.fna"
                  onChange={handleFileChange}
                />
                <label
                  className={`
                    group flex min-h-14 cursor-pointer items-center gap-3 rounded-[8px] border bg-white px-3
                    shadow-sm shadow-stone-900/5 transition hover:border-emerald-500 hover:bg-emerald-50/45
                    ${file ? "border-emerald-300" : "border-dashed border-emerald-200"}
                  `}
                  htmlFor="fasta-file"
                >
                  <span
                    className={`
                      flex size-9 shrink-0 items-center justify-center rounded-[7px] border transition
                      ${file ? "border-emerald-200 bg-emerald-50 text-emerald-700" : "border-emerald-100 bg-emerald-50/50 text-emerald-700"}
                    `}
                  >
                    <UploadCloud className="size-4"/>
                  </span>
                  <span className="min-w-0 flex-1">
                    <span className="block truncate text-sm font-medium text-stone-900">
                      {file ? file.name : "Attach FASTA file"}
                    </span>
                    <span className="block text-xs text-stone-500">
                      {file ? `${formatFileSize(file.size)} selected` : "Accepted formats: .fasta, .fa, .fna"}
                    </span>
                  </span>
                </label>
              </div>

            </form>
          </section>
        ) : null}

        {isLoading ? (
          <section className="flex min-h-[78vh] items-center justify-center">
            <div className="w-full max-w-lg">
              <div className="flex items-start justify-between gap-4">
                <div className="flex min-w-0 items-center gap-3">
                  <span className="flex size-11 shrink-0 items-center justify-center rounded-[8px] border border-emerald-100 bg-white text-emerald-700 shadow-sm shadow-emerald-900/5">
                    <Loader2 className="size-5 animate-spin"/>
                  </span>
                  <div className="min-w-0">
                    <p className="text-sm font-medium text-stone-950">Scanning motif</p>
                    <p className="mt-1 truncate text-sm text-stone-500">{loadingStage}</p>
                  </div>
                </div>
                <span className="font-mono text-3xl font-semibold tabular-nums text-emerald-800">
                  {loadingProgress}%
                </span>
              </div>

              <div className="mt-6 h-2 overflow-hidden rounded-full bg-emerald-100">
                <div
                  className="h-2 rounded-full bg-emerald-700 transition-[width] duration-300 ease-out"
                  style={{width: `${loadingProgress}%`}}
                />
              </div>

              <div className="mt-5 rounded-[8px] border border-emerald-100 bg-white p-3 shadow-sm shadow-emerald-900/5">
                <div className="flex items-center justify-between gap-3 text-xs text-stone-500">
                  <span>Motif</span>
                  <span className="font-mono font-medium text-emerald-800">
                    {motifBaseCount} {motifBaseCount === 1 ? "base" : "bases"}
                  </span>
                </div>
                <p className="mt-2 max-h-16 overflow-y-auto break-all font-mono text-sm font-medium leading-5 text-stone-950">
                  {motifPattern}
                </p>
              </div>

              <div className="mt-3 flex min-h-14 items-center gap-3 rounded-[8px] border border-emerald-100 bg-white px-3 shadow-sm shadow-emerald-900/5">
                <span className="flex size-9 shrink-0 items-center justify-center rounded-[7px] bg-emerald-50 text-emerald-700">
                  <UploadCloud className="size-4"/>
                </span>
                <span className="min-w-0 flex-1">
                  <span className="block truncate text-sm font-medium text-stone-900">
                    {file?.name ?? "FASTA file"}
                  </span>
                  <span className="block text-xs text-stone-500">
                    {profile
                      ? `${formatInteger(profile.records)} records, ${formatInteger(profile.canonicalBases)} canonical bases`
                      : "Processing FASTA"}
                  </span>
                </span>
              </div>
            </div>
          </section>
        ) : null}

        {result ? (
          <section className="rounded-[10px] border border-stone-200 bg-white shadow-sm">
            <div className="border-b border-stone-200 p-4">
              <div className="mt-1 flex flex-wrap items-center justify-between gap-3">
                <h2 className="text-xl font-semibold">Results</h2>
                <div className="flex gap-2">
                  <Button className="bg-white" type="button" variant="outline" onClick={() => setResult(null)}>
                    New scan
                  </Button>
                  <Button className="bg-white" type="button" variant="outline" onClick={exportCsv}>
                    <FileSpreadsheet className="size-4"/>
                    CSV
                  </Button>
                  <Button className="bg-white" type="button" variant="outline" onClick={exportJson}>
                    <FileJson className="size-4"/>
                    JSON
                  </Button>
                </div>
              </div>
            </div>

            <div className="p-4">
              <div className="grid gap-3 sm:grid-cols-2 lg:grid-cols-4">
                <MetricCard label="Total matches" value={formatInteger(stats.totalMatches)}/>
                <MetricCard label="Genes with hits" value={formatInteger(stats.uniqueGenes)}/>
                <MetricCard label="Best p-value"
                            value={stats.bestPValue === null ? "-" : formatPValue(stats.bestPValue)}/>
                <MetricCard label="Max score" value={stats.maxScore === null ? "-" : String(stats.maxScore)}/>
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
              </div>

              {activeTab === "overview" ? (
                <OverviewPanel metadata={result.metadata} profile={profile} stats={stats}
                               completedAt={result.completedAt}/>
              ) : null}
              {activeTab === "genes" ? <GenesPanel rows={geneSummary}/> : null}
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
            </div>
          </section>
        ) : null}
      </div>
    </main>
  )
}

function MotifPreview({columns, motifPattern}: { columns: Base[][]; motifPattern: string }) {
  if (!motifPattern.trim() || !columns.length) {
    return null
  }

  return (
    <div>
      <div className="overflow-x-auto rounded-[8px] border border-emerald-100 bg-white p-2 shadow-sm shadow-emerald-900/5">
        <div className="flex w-max min-w-full gap-1.5">
          {columns.map((column, idx) => (
            <div
              className="group flex min-w-9 flex-col items-center gap-1 rounded-[6px] px-1.5 py-1.5 transition hover:bg-stone-50"
              key={`${motifPattern}-${idx}`}
              title={`Position ${idx + 1}`}
            >
              <span className="font-mono text-[10px] text-stone-400">{idx + 1}</span>
              <div className="flex min-h-5 flex-col gap-1">
                {column.map((base) => (
                  <span
                    className={`flex size-5 items-center justify-center rounded-[5px] border text-[11px] font-bold leading-none ${BASE_SWATCHES[base]}`}
                    key={`${idx}-${base}`}
                  >
                    {base}
                  </span>
                ))}
              </div>
            </div>
          ))}
        </div>
      </div>
    </div>
  )
}

function MetricCard({label, value}: { label: string; value: string }) {
  return (
    <article className="rounded-[8px] border border-stone-200 bg-white p-3">
      <p className="text-xs font-medium uppercase text-stone-600">{label}</p>
      <p className="mt-1 font-mono text-lg font-semibold text-stone-950">{value}</p>
    </article>
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
          <StatLine label="Input file" value={metadata.filename}/>
          <StatLine label="Motif pattern" value={metadata.motif_pattern}/>
          <StatLine label="Motif length" value={String(metadata.motif_length)}/>
          <StatLine label="Cutoff mode" value={metadata.cutoff_mode}/>
          <StatLine
            label="Score threshold"
            value={metadata.score_threshold === null ? "derived from p-value" : String(metadata.score_threshold)}
          />
          <StatLine label="P-value threshold" value={DEFAULT_P_VALUE_THRESHOLD}/>
          <StatLine label="Completed at" value={new Date(completedAt).toLocaleString()}/>
          <StatLine label="Forward source hits" value={formatInteger(stats.forwardSource)}/>
          <StatLine label="Reverse source hits" value={formatInteger(stats.reverseSource)}/>
          <StatLine label="Canonical bases in file" value={profile ? formatInteger(profile.canonicalBases) : "-"}/>
          <StatLine
            label="P-value method"
            value={metadata.pvalue_engine?.method ?? "discretized distribution"}
          />
          <StatLine
            label="Granularity levels"
            value={String(metadata.pvalue_engine?.levels_used ?? "-")}
          />
        </dl>
      </div>

      <div className="rounded-[8px] border border-stone-200 bg-stone-50 p-4">
        <h3 className="text-sm font-semibold text-stone-900">Background frequencies</h3>
        <div className="mt-3 space-y-3">
          {BASES.map((base) => (
            <div className="grid grid-cols-[24px_minmax(0,1fr)_50px] items-center gap-2" key={base}>
              <span
                className={`flex size-6 items-center justify-center rounded-[6px] border text-xs font-bold ${BASE_SWATCHES[base]}`}
              >
                {base}
              </span>
              <div className="h-2 rounded-full bg-stone-200">
                <div
                  className="h-2 rounded-full bg-emerald-700"
                  style={{width: `${Math.max(2, (metadata.background_frequencies[base] ?? 0) * 100)}%`}}
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

function GenesPanel({rows}: { rows: GeneSummaryRow[] }) {
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
          <Search className="pointer-events-none absolute left-3 top-1/2 size-4 -translate-y-1/2 text-stone-500"/>
          <input
            className={`${CONTROL_CLASS} pl-9`}
            placeholder="Search gene or sequence"
            value={searchTerm}
            onChange={(event) => setSearchTerm(event.target.value)}
          />
        </label>
        <label className="relative block">
          <Filter className="pointer-events-none absolute left-3 top-1/2 size-4 -translate-y-1/2 text-stone-500"/>
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
          <ArrowDownUp className="pointer-events-none absolute left-3 top-1/2 size-4 -translate-y-1/2 text-stone-500"/>
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
        <table className="min-w-[960px] w-full text-left text-sm">
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
                <SequenceSwatches sequence={match.matched_sequence}/>
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

function StatLine({label, value}: { label: string; value: string }) {
  return (
    <div className="flex items-center justify-between gap-3">
      <dt className="truncate text-stone-500">{label}</dt>
      <dd className="truncate text-right font-mono text-stone-900">{value}</dd>
    </div>
  )
}

function SequenceSwatches({sequence}: { sequence: string }) {
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
  const counts: Record<Base, number> = {A: 0, C: 0, G: 0, T: 0}
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

async function readScanStream(
  response: Response,
  onProgress: (progress: number, stage: string) => void,
): Promise<SingleEnginePayload> {
  if (!response.body) {
    return parseScanResponse(await response.text())
  }

  const reader = response.body.getReader()
  const decoder = new TextDecoder()
  let buffer = ""
  let payload: SingleEnginePayload | null = null

  const consumeLine = (line: string) => {
    const trimmed = line.trim()
    if (!trimmed) return

    const event = JSON.parse(trimmed) as ScanStreamEvent
    if (event.type === "progress") {
      onProgress(normalizeProgress(event.progress), event.stage)
      return
    }

    if (event.type === "error") {
      throw new Error(event.message)
    }

    if (event.type === "result") {
      payload = event.payload
    }
  }

  while (true) {
    const {value, done} = await reader.read()
    buffer += decoder.decode(value ?? new Uint8Array(), {stream: !done})

    const lines = buffer.split(/\r?\n/)
    buffer = lines.pop() ?? ""
    for (const line of lines) consumeLine(line)

    if (done) break
  }

  if (buffer.trim()) consumeLine(buffer)
  if (!payload) throw new Error("API stream finished without scan results.")
  return payload
}

function normalizeProgress(progress: number) {
  if (!Number.isFinite(progress)) return 0
  return Math.max(0, Math.min(100, Math.round(progress)))
}

function parseScanResponse(raw: string): SingleEnginePayload {
  const trimmed = raw.trim()
  if (!trimmed) throw new Error("API returned an empty response.")

  try {
    return JSON.parse(trimmed) as SingleEnginePayload
  } catch {
    const payloads = trimmed.split(/\r?\n/).filter(Boolean).map((line) => JSON.parse(line) as SingleEnginePayload)
    const payload = payloads.find((candidate) => Boolean(candidate.metadata && candidate.matches))
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

function formatFileSize(bytes: number) {
  if (!Number.isFinite(bytes) || bytes <= 0) return "0 B"
  const units = ["B", "KB", "MB", "GB"]
  const exponent = Math.min(Math.floor(Math.log(bytes) / Math.log(1024)), units.length - 1)
  const value = bytes / 1024 ** exponent
  return `${value.toFixed(value >= 10 || exponent === 0 ? 0 : 1)} ${units[exponent]}`
}

function safeTimestamp() {
  return new Date().toISOString().replaceAll(":", "-").replace(/\.\d{3}Z$/, "Z")
}

function escapeCsv(value: string) {
  return `"${value.replaceAll('"', '""')}"`
}

function downloadText(filename: string, text: string, type: string) {
  const blob = new Blob([text], {type})
  const url = URL.createObjectURL(blob)
  const anchor = document.createElement("a")
  anchor.href = url
  anchor.download = filename
  anchor.click()
  URL.revokeObjectURL(url)
}

export default App
