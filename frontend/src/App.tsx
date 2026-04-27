import {
  Activity,
  ArrowDownUp,
  BarChart3,
  CheckCircle2,
  ClipboardList,
  Dna,
  Download,
  FileJson,
  FileSpreadsheet,
  FlaskConical,
  Loader2,
  Microscope,
  Play,
  Search,
  Settings2,
  SlidersHorizontal,
  Table2,
  TriangleAlert,
  UploadCloud,
} from "lucide-react"
import { useMemo, useState } from "react"
import type { ChangeEvent, FormEvent, ReactNode } from "react"

import { Button } from "@/components/ui/button"

type Base = "A" | "C" | "G" | "T"
type RunMode = "pvalue" | "score"
type ActiveView = "table" | "summary" | "audit"
type StrandFilter = "all" | "forward" | "reverse"
type SortMode = "pvalue" | "score" | "gene" | "position"

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
  bestHit: MotifMatch | null
  maxScore: number | null
  medianPValue: number | null
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

const BASE_HEX: Record<Base, string> = {
  A: "#059669",
  C: "#0284c7",
  G: "#d97706",
  T: "#e11d48",
}

const CONTROL_CLASS =
  "h-10 w-full rounded-[8px] border border-stone-300 bg-white px-3 text-sm text-stone-950 shadow-sm outline-none transition focus:border-emerald-600 focus:ring-3 focus:ring-emerald-600/15"

const TEXTAREA_CLASS =
  "min-h-20 w-full resize-y rounded-[8px] border border-stone-300 bg-white px-3 py-2 text-sm text-stone-950 shadow-sm outline-none transition focus:border-emerald-600 focus:ring-3 focus:ring-emerald-600/15"

function App() {
  const [file, setFile] = useState<File | null>(null)
  const [profile, setProfile] = useState<FastaProfile | null>(null)
  const [hypothesis, setHypothesis] = useState("")
  const [motifPattern, setMotifPattern] = useState("TATAWA")
  const [runMode, setRunMode] = useState<RunMode>("pvalue")
  const [pValueThreshold, setPValueThreshold] = useState("1e-4")
  const [scoreThreshold, setScoreThreshold] = useState("")
  const [scanReverseComplement, setScanReverseComplement] = useState(true)
  const [isLoading, setIsLoading] = useState(false)
  const [error, setError] = useState("")
  const [result, setResult] = useState<ScanResult | null>(null)
  const [activeView, setActiveView] = useState<ActiveView>("table")
  const [searchTerm, setSearchTerm] = useState("")
  const [strandFilter, setStrandFilter] = useState<StrandFilter>("all")
  const [sortMode, setSortMode] = useState<SortMode>("pvalue")

  const motifColumns = useMemo(() => parseMotifPattern(motifPattern), [motifPattern])
  const matches = useMemo(() => result?.matches ?? [], [result])

  const filteredMatches = useMemo(() => {
    const query = searchTerm.trim().toLowerCase()
    const filtered = matches.filter((match) => {
      const matchesQuery =
        !query ||
        match.gene.toLowerCase().includes(query) ||
        match.matched_sequence.toLowerCase().includes(query)
      const matchesStrand = strandFilter === "all" || match.source_strand === strandFilter

      return matchesQuery && matchesStrand
    })

    return [...filtered].sort((a, b) => sortMatches(a, b, sortMode))
  }, [matches, searchTerm, sortMode, strandFilter])

  const visibleMatches = filteredMatches.slice(0, 300)
  const stats = useMemo(() => calculateStats(matches), [matches])
  const geneSummary = useMemo(() => summarizeGenes(matches), [matches])

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
      setError("Unable to read the FASTA file in the browser.")
    }
  }

  async function handleSubmit(event: FormEvent<HTMLFormElement>) {
    event.preventDefault()
    setError("")

    if (!file) {
      setError("Select a FASTA file before running the scan.")
      return
    }

    const trimmedMotif = motifPattern.trim()
    if (!trimmedMotif) {
      setError("Define a motif pattern before running the scan.")
      return
    }

    if (runMode === "pvalue") {
      const pValue = Number(pValueThreshold)
      if (!Number.isFinite(pValue) || pValue <= 0 || pValue > 1) {
        setError("P-value threshold must be greater than 0 and less than or equal to 1.")
        return
      }
    }

    if (runMode === "score") {
      const score = Number(scoreThreshold)
      if (!Number.isInteger(score)) {
        setError("Score threshold must be an integer because the backend uses an integer PWM.")
        return
      }
    }

    const formData = new FormData()
    formData.append("file", file)
    formData.append("motif_pattern", trimmedMotif)
    formData.append("scan_reverse_complement", String(scanReverseComplement))

    if (runMode === "pvalue") {
      formData.append("p_value_threshold", pValueThreshold)
    } else {
      formData.append("score_threshold", scoreThreshold)
    }

    setIsLoading(true)

    try {
      const response = await fetch(`${API_BASE}/motifs/`, {
        method: "POST",
        body: formData,
      })
      const raw = await response.text()

      if (!response.ok) {
        throw new Error(readApiError(raw))
      }

      const payload = parseScanResponse(raw)
      setResult({
        metadata: payload.metadata,
        matches: payload.matches.map((match) => ({
          ...match,
          negLog10PValue: match["-log10_pval"],
        })),
        completedAt: new Date().toISOString(),
      })
      setActiveView("table")
    } catch (caughtError) {
      setResult(null)
      setError(caughtError instanceof Error ? caughtError.message : "The motif scan failed.")
    } finally {
      setIsLoading(false)
    }
  }

  function exportJson() {
    if (!result) return

    downloadText(
      `motif-scan-${safeTimestamp()}.json`,
      JSON.stringify(
        {
          hypothesis,
          fasta_profile: profile,
          ...result,
        },
        null,
        2,
      ),
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
      "neg_log10_p_value",
      "record_orientation",
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
      match.record_orientation,
    ])

    downloadText(
      `motif-matches-${safeTimestamp()}.csv`,
      [header, ...rows].map((row) => row.map(escapeCsv).join(",")).join("\n"),
      "text/csv",
    )
  }

  return (
    <main className="min-h-screen bg-[#f5f8f4] text-stone-950">
      <header className="border-b border-stone-200 bg-white/90">
        <div className="mx-auto flex max-w-[1440px] flex-wrap items-center justify-between gap-4 px-4 py-4 sm:px-6 lg:px-8">
          <div className="flex items-center gap-3">
            <div className="flex size-10 items-center justify-center rounded-[8px] border border-emerald-200 bg-emerald-50 text-emerald-800">
              <Microscope className="size-5" />
            </div>
            <div>
              <p className="text-xs font-semibold uppercase text-emerald-700">Motif evidence workbench</p>
              <h1 className="text-xl font-semibold text-stone-950 sm:text-2xl">
                Position weight matrix screening
              </h1>
            </div>
          </div>
          <div className="flex max-w-full items-center gap-2 rounded-[8px] border border-stone-200 bg-stone-50 px-3 py-2 text-xs text-stone-600">
            <Activity className="size-4 text-emerald-700" />
            <span className="shrink-0 font-medium text-stone-800">API</span>
            <code className="truncate text-stone-600">{API_BASE}/motifs/</code>
          </div>
        </div>
      </header>

      <div className="mx-auto grid max-w-[1440px] gap-5 px-4 py-5 sm:px-6 lg:grid-cols-[390px_minmax(0,1fr)] lg:px-8">
        <form className="space-y-5" onSubmit={handleSubmit}>
          <section className="rounded-[8px] border border-stone-200 bg-white p-4 shadow-sm">
            <SectionHeading
              icon={<FlaskConical className="size-4" />}
              kicker="Study frame"
              title="Hypothesis and FASTA dataset"
            />

            <label className="mt-4 block text-sm font-medium text-stone-800" htmlFor="hypothesis">
              Working hypothesis
            </label>
            <textarea
              id="hypothesis"
              className={`${TEXTAREA_CLASS} mt-2`}
              value={hypothesis}
              onChange={(event) => setHypothesis(event.target.value)}
              placeholder="Example: Candidate promoter motif frequency differs from the uploaded genomic background."
            />

            <div className="mt-4">
              <input
                className="sr-only"
                id="fasta-file"
                type="file"
                accept=".fasta,.fa,.fna"
                onChange={handleFileChange}
              />
              <label
                className="flex min-h-28 cursor-pointer flex-col items-center justify-center gap-2 rounded-[8px] border border-dashed border-emerald-300 bg-emerald-50/70 px-4 py-5 text-center transition hover:border-emerald-500 hover:bg-emerald-50"
                htmlFor="fasta-file"
              >
                <UploadCloud className="size-6 text-emerald-700" />
                <span className="max-w-full truncate text-sm font-semibold text-stone-900">
                  {file ? file.name : "Choose FASTA file"}
                </span>
                <span className="text-xs text-stone-600">.fasta, .fa, or .fna</span>
              </label>
            </div>

            {profile ? <FastaProfilePanel profile={profile} /> : null}
          </section>

          <section className="rounded-[8px] border border-stone-200 bg-white p-4 shadow-sm">
            <SectionHeading
              icon={<Dna className="size-4" />}
              kicker="Motif model"
              title="Pattern and decision rule"
            />

            <label className="mt-4 block text-sm font-medium text-stone-800" htmlFor="motif-pattern">
              IUPAC motif pattern
            </label>
            <input
              id="motif-pattern"
              className={`${CONTROL_CLASS} mt-2 font-mono uppercase`}
              value={motifPattern}
              onChange={(event) => setMotifPattern(event.target.value.toUpperCase())}
              placeholder="TATAWA"
            />
            <MotifPreview columns={motifColumns} />

            <div className="mt-5">
              <span className="block text-sm font-medium text-stone-800">Threshold mode</span>
              <div className="mt-2 grid grid-cols-2 gap-1 rounded-[8px] border border-stone-200 bg-stone-100 p-1">
                <button
                  type="button"
                  className={segmentClass(runMode === "pvalue")}
                  onClick={() => setRunMode("pvalue")}
                >
                  P-value
                </button>
                <button
                  type="button"
                  className={segmentClass(runMode === "score")}
                  onClick={() => setRunMode("score")}
                >
                  Score
                </button>
              </div>
            </div>

            {runMode === "pvalue" ? (
              <div className="mt-4">
                <label className="block text-sm font-medium text-stone-800" htmlFor="p-value">
                  P-value threshold
                </label>
                <input
                  id="p-value"
                  className={`${CONTROL_CLASS} mt-2 font-mono`}
                  value={pValueThreshold}
                  onChange={(event) => setPValueThreshold(event.target.value)}
                  inputMode="decimal"
                />
              </div>
            ) : (
              <div className="mt-4">
                <label className="block text-sm font-medium text-stone-800" htmlFor="score-threshold">
                  Integer score threshold
                </label>
                <input
                  id="score-threshold"
                  className={`${CONTROL_CLASS} mt-2 font-mono`}
                  value={scoreThreshold}
                  onChange={(event) => setScoreThreshold(event.target.value)}
                  inputMode="numeric"
                  placeholder="120"
                />
              </div>
            )}

            <label className="mt-4 flex items-center justify-between gap-3 rounded-[8px] border border-stone-200 bg-stone-50 px-3 py-2">
              <span className="text-sm font-medium text-stone-800">Scan reverse complement</span>
              <input
                className="size-5 accent-emerald-700"
                type="checkbox"
                checked={scanReverseComplement}
                onChange={(event) => setScanReverseComplement(event.target.checked)}
              />
            </label>

            <Button className="mt-5 h-10 w-full bg-emerald-700 text-white hover:bg-emerald-800" disabled={isLoading}>
              {isLoading ? <Loader2 className="size-4 animate-spin" /> : <Play className="size-4" />}
              {isLoading ? "Scanning" : "Run motif scan"}
            </Button>
          </section>

          {error ? (
            <div className="flex gap-3 rounded-[8px] border border-rose-200 bg-rose-50 p-3 text-sm text-rose-900">
              <TriangleAlert className="mt-0.5 size-4 shrink-0" />
              <span>{error}</span>
            </div>
          ) : null}
        </form>

        <section className="min-w-0 space-y-5">
          <div className="grid gap-3 sm:grid-cols-2 xl:grid-cols-4">
            <MetricCard
              icon={<ClipboardList className="size-4" />}
              label="Matches"
              value={formatInteger(stats.totalMatches)}
              detail={result ? "Current result set" : "Awaiting scan"}
            />
            <MetricCard
              icon={<Dna className="size-4" />}
              label="Genes with evidence"
              value={formatInteger(stats.uniqueGenes)}
              detail={result ? "Unique FASTA identifiers" : "No result yet"}
            />
            <MetricCard
              icon={<ArrowDownUp className="size-4" />}
              label="Source strands"
              value={`${stats.forwardSource}/${stats.reverseSource}`}
              detail="Forward / reverse"
            />
            <MetricCard
              icon={<CheckCircle2 className="size-4" />}
              label="Best p-value"
              value={stats.bestHit ? formatPValue(stats.bestHit.p_value_numeric) : "—"}
              detail={stats.bestHit?.gene ?? "No hit selected"}
            />
          </div>

          <section className="rounded-[8px] border border-stone-200 bg-white p-4 shadow-sm">
            <div className="flex flex-wrap items-center justify-between gap-3">
              <SectionHeading
                icon={<BarChart3 className="size-4" />}
                kicker="Evidence landscape"
                title="Sorted motif-hit significance"
              />
              <div className="flex gap-2">
                <Button
                  type="button"
                  variant="outline"
                  className="bg-white"
                  disabled={!result}
                  onClick={exportCsv}
                >
                  <FileSpreadsheet className="size-4" />
                  CSV
                </Button>
                <Button
                  type="button"
                  variant="outline"
                  className="bg-white"
                  disabled={!result}
                  onClick={exportJson}
                >
                  <FileJson className="size-4" />
                  JSON
                </Button>
              </div>
            </div>
            <EvidenceChart matches={filteredMatches} />
          </section>

          <section className="rounded-[8px] border border-stone-200 bg-white shadow-sm">
            <div className="flex flex-wrap items-center justify-between gap-3 border-b border-stone-200 p-4">
              <div className="flex flex-wrap gap-2">
                <ViewButton active={activeView === "table"} onClick={() => setActiveView("table")}>
                  <Table2 className="size-4" />
                  Evidence table
                </ViewButton>
                <ViewButton active={activeView === "summary"} onClick={() => setActiveView("summary")}>
                  <BarChart3 className="size-4" />
                  Locus summary
                </ViewButton>
                <ViewButton active={activeView === "audit"} onClick={() => setActiveView("audit")}>
                  <Settings2 className="size-4" />
                  Method audit
                </ViewButton>
              </div>
              <Button
                type="button"
                variant="outline"
                className="bg-white"
                disabled={!result}
                onClick={exportJson}
              >
                <Download className="size-4" />
                Export run
              </Button>
            </div>

            {activeView === "table" ? (
              <EvidenceTable
                filteredMatches={filteredMatches}
                visibleMatches={visibleMatches}
                searchTerm={searchTerm}
                setSearchTerm={setSearchTerm}
                sortMode={sortMode}
                setSortMode={setSortMode}
                strandFilter={strandFilter}
                setStrandFilter={setStrandFilter}
              />
            ) : null}

            {activeView === "summary" ? (
              <LocusSummary rows={geneSummary} stats={stats} />
            ) : null}

            {activeView === "audit" ? (
              <MethodAudit
                hypothesis={hypothesis}
                metadata={result?.metadata ?? null}
                profile={profile}
                runMode={runMode}
                pValueThreshold={pValueThreshold}
                scoreThreshold={scoreThreshold}
                scanReverseComplement={scanReverseComplement}
                completedAt={result?.completedAt ?? ""}
              />
            ) : null}
          </section>
        </section>
      </div>
    </main>
  )
}

function SectionHeading({
  icon,
  kicker,
  title,
}: {
  icon: ReactNode
  kicker: string
  title: string
}) {
  return (
    <div className="flex items-start gap-3">
      <div className="mt-0.5 flex size-8 items-center justify-center rounded-[8px] border border-stone-200 bg-stone-50 text-emerald-700">
        {icon}
      </div>
      <div>
        <p className="text-xs font-semibold uppercase text-emerald-700">{kicker}</p>
        <h2 className="text-base font-semibold text-stone-950">{title}</h2>
      </div>
    </div>
  )
}

function FastaProfilePanel({ profile }: { profile: FastaProfile }) {
  return (
    <div className="mt-4 rounded-[8px] border border-stone-200 bg-stone-50 p-3">
      <div className="flex items-center justify-between gap-3 text-sm">
        <span className="font-medium text-stone-800">Dataset profile</span>
        <span className="text-xs text-stone-600">{formatBytes(profile.sizeBytes)}</span>
      </div>
      <dl className="mt-3 grid grid-cols-2 gap-2 text-sm">
        <ProfileItem label="Records" value={formatInteger(profile.records)} />
        <ProfileItem label="Bases" value={formatInteger(profile.totalBases)} />
        <ProfileItem label="GC" value={formatPercent(profile.gcRatio)} />
        <ProfileItem label="Ambiguous" value={formatInteger(profile.ambiguousBases)} />
      </dl>
      {profile.firstHeader ? (
        <p className="mt-3 truncate border-t border-stone-200 pt-3 text-xs text-stone-600">{profile.firstHeader}</p>
      ) : null}
    </div>
  )
}

function ProfileItem({ label, value }: { label: string; value: string }) {
  return (
    <div className="rounded-[8px] border border-stone-200 bg-white px-3 py-2">
      <dt className="text-xs text-stone-500">{label}</dt>
      <dd className="mt-1 font-semibold text-stone-950">{value}</dd>
    </div>
  )
}

function MotifPreview({ columns }: { columns: Base[][] }) {
  if (!columns.length) {
    return (
      <div className="mt-3 rounded-[8px] border border-amber-200 bg-amber-50 px-3 py-2 text-sm text-amber-900">
        Motif preview is waiting for a valid IUPAC pattern.
      </div>
    )
  }

  return (
    <div className="mt-3 overflow-x-auto rounded-[8px] border border-stone-200 bg-stone-50 p-2">
      <div className="flex min-w-max gap-1">
        {columns.map((column, index) => (
          <div
            className="flex w-11 flex-col items-center gap-1 rounded-[8px] border border-stone-200 bg-white p-1"
            key={`${column.join("")}-${index}`}
            title={`Position ${index + 1}`}
          >
            <span className="text-[10px] font-medium text-stone-500">{index + 1}</span>
            <div className="grid grid-cols-1 gap-1">
              {column.map((base) => (
                <span
                  className={`flex size-6 items-center justify-center rounded-[6px] border text-xs font-bold ${BASE_SWATCHES[base]}`}
                  key={base}
                >
                  {base}
                </span>
              ))}
            </div>
          </div>
        ))}
      </div>
    </div>
  )
}

function MetricCard({
  icon,
  label,
  value,
  detail,
}: {
  icon: ReactNode
  label: string
  value: string
  detail: string
}) {
  return (
    <div className="rounded-[8px] border border-stone-200 bg-white p-4 shadow-sm">
      <div className="flex items-center justify-between gap-3">
        <span className="text-sm font-medium text-stone-600">{label}</span>
        <span className="flex size-8 items-center justify-center rounded-[8px] bg-emerald-50 text-emerald-700">
          {icon}
        </span>
      </div>
      <p className="mt-3 truncate text-2xl font-semibold text-stone-950">{value}</p>
      <p className="mt-1 truncate text-xs text-stone-500">{detail}</p>
    </div>
  )
}

function EvidenceChart({ matches }: { matches: MotifMatch[] }) {
  const topMatches = matches
    .filter((match) => Number.isFinite(match.p_value_numeric) && match.p_value_numeric > 0)
    .slice(0, 48)
  const values = topMatches.map((match) => match.negLog10PValue ?? -Math.log10(match.p_value_numeric))
  const maxValue = Math.max(1, ...values)
  const chartWidth = 760
  const chartHeight = 210
  const plotTop = 18
  const plotBottom = 172
  const barGap = 3
  const barWidth = topMatches.length ? Math.max(4, (chartWidth - 64) / topMatches.length - barGap) : 0

  if (!topMatches.length) {
    return (
      <div className="mt-4 flex min-h-48 items-center justify-center rounded-[8px] border border-dashed border-stone-300 bg-stone-50 text-sm text-stone-600">
        No scored motif evidence to plot yet.
      </div>
    )
  }

  return (
    <div className="mt-4 overflow-x-auto rounded-[8px] border border-stone-200 bg-stone-50 p-2">
      <svg aria-label="Sorted motif hit significance chart" className="h-56 min-w-[760px] w-full" viewBox={`0 0 ${chartWidth} ${chartHeight}`}>
        <line x1="42" x2="736" y1={plotBottom} y2={plotBottom} stroke="#d6d3d1" strokeWidth="1" />
        <line x1="42" x2="42" y1={plotTop} y2={plotBottom} stroke="#d6d3d1" strokeWidth="1" />
        <text fill="#78716c" fontSize="11" x="8" y="22">
          -log10 p
        </text>
        <text fill="#78716c" fontSize="11" x="650" y="197">
          top hits by significance
        </text>
        {topMatches.map((match, index) => {
          const value = values[index]
          const height = Math.max(2, (value / maxValue) * (plotBottom - plotTop))
          const x = 48 + index * (barWidth + barGap)
          const y = plotBottom - height
          const fill = match.source_strand === "reverse" ? "#be123c" : "#047857"

          return (
            <g key={`${match.gene}-${match.source_position}-${index}`}>
              <rect fill={fill} height={height} rx="2" width={barWidth} x={x} y={y} />
              {index % 8 === 0 ? (
                <text fill="#78716c" fontSize="9" textAnchor="middle" x={x + barWidth / 2} y="187">
                  {index + 1}
                </text>
              ) : null}
            </g>
          )
        })}
      </svg>
    </div>
  )
}

function ViewButton({
  active,
  children,
  onClick,
}: {
  active: boolean
  children: ReactNode
  onClick: () => void
}) {
  return (
    <button
      className={`inline-flex h-9 items-center gap-2 rounded-[8px] border px-3 text-sm font-medium transition ${
        active
          ? "border-emerald-700 bg-emerald-700 text-white"
          : "border-stone-200 bg-white text-stone-700 hover:border-emerald-300 hover:text-emerald-800"
      }`}
      type="button"
      onClick={onClick}
    >
      {children}
    </button>
  )
}

function EvidenceTable({
  filteredMatches,
  visibleMatches,
  searchTerm,
  setSearchTerm,
  sortMode,
  setSortMode,
  strandFilter,
  setStrandFilter,
}: {
  filteredMatches: MotifMatch[]
  visibleMatches: MotifMatch[]
  searchTerm: string
  setSearchTerm: (value: string) => void
  sortMode: SortMode
  setSortMode: (value: SortMode) => void
  strandFilter: StrandFilter
  setStrandFilter: (value: StrandFilter) => void
}) {
  return (
    <div className="p-4">
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
          <SlidersHorizontal className="pointer-events-none absolute left-3 top-1/2 size-4 -translate-y-1/2 text-stone-500" />
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

      <div className="mt-4 overflow-x-auto rounded-[8px] border border-stone-200">
        <table className="min-w-[920px] w-full border-collapse text-left text-sm">
          <thead className="bg-stone-100 text-xs uppercase text-stone-600">
            <tr>
              <th className="px-3 py-3 font-semibold">Gene</th>
              <th className="px-3 py-3 font-semibold">Matched sequence</th>
              <th className="px-3 py-3 font-semibold">Coordinates</th>
              <th className="px-3 py-3 font-semibold">Strand</th>
              <th className="px-3 py-3 font-semibold">Score</th>
              <th className="px-3 py-3 font-semibold">P-value</th>
              <th className="px-3 py-3 font-semibold">-log10 p</th>
            </tr>
          </thead>
          <tbody className="divide-y divide-stone-200 bg-white">
            {visibleMatches.map((match, index) => (
              <tr className="hover:bg-emerald-50/50" key={`${match.gene}-${match.source_position}-${index}`}>
                <td className="max-w-52 px-3 py-3">
                  <span className="block truncate font-medium text-stone-950" title={match.gene}>
                    {match.gene}
                  </span>
                  <span className="text-xs text-stone-500">{match.record_orientation}</span>
                </td>
                <td className="px-3 py-3">
                  <SequenceSwatches sequence={match.matched_sequence} />
                </td>
                <td className="px-3 py-3 font-mono text-xs text-stone-700">
                  {match.source_position}–{match.source_end_position}
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
                <td className="px-3 py-3 font-mono text-stone-800">{match.score}</td>
                <td className="px-3 py-3 font-mono text-stone-800">{match.p_value}</td>
                <td className="px-3 py-3 font-mono text-stone-800">
                  {match.negLog10PValue?.toFixed(2) ?? "—"}
                </td>
              </tr>
            ))}
            {!visibleMatches.length ? (
              <tr>
                <td className="px-3 py-12 text-center text-sm text-stone-500" colSpan={7}>
                  No motif hits match the current view.
                </td>
              </tr>
            ) : null}
          </tbody>
        </table>
      </div>

      <p className="mt-3 text-xs text-stone-500">
        Showing {formatInteger(visibleMatches.length)} of {formatInteger(filteredMatches.length)} filtered hits.
      </p>
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

function LocusSummary({ rows, stats }: { rows: GeneSummaryRow[]; stats: RunStats }) {
  return (
    <div className="p-4">
      <div className="grid gap-3 md:grid-cols-3">
        <SummaryBlock label="Median p-value" value={stats.medianPValue ? formatPValue(stats.medianPValue) : "—"} />
        <SummaryBlock label="Maximum score" value={stats.maxScore === null ? "—" : String(stats.maxScore)} />
        <SummaryBlock label="Gene count" value={formatInteger(rows.length)} />
      </div>

      <div className="mt-4 overflow-x-auto rounded-[8px] border border-stone-200">
        <table className="min-w-[760px] w-full text-left text-sm">
          <thead className="bg-stone-100 text-xs uppercase text-stone-600">
            <tr>
              <th className="px-3 py-3">Gene</th>
              <th className="px-3 py-3">Hits</th>
              <th className="px-3 py-3">Best p-value</th>
              <th className="px-3 py-3">Best score</th>
              <th className="px-3 py-3">First position</th>
              <th className="px-3 py-3">Source strands</th>
            </tr>
          </thead>
          <tbody className="divide-y divide-stone-200 bg-white">
            {rows.slice(0, 200).map((row) => (
              <tr className="hover:bg-emerald-50/50" key={row.gene}>
                <td className="max-w-72 truncate px-3 py-3 font-medium text-stone-950">{row.gene}</td>
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
                  No loci to summarize yet.
                </td>
              </tr>
            ) : null}
          </tbody>
        </table>
      </div>
    </div>
  )
}

function SummaryBlock({ label, value }: { label: string; value: string }) {
  return (
    <div className="rounded-[8px] border border-stone-200 bg-stone-50 p-3">
      <p className="text-xs font-medium uppercase text-stone-500">{label}</p>
      <p className="mt-2 text-xl font-semibold text-stone-950">{value}</p>
    </div>
  )
}

function MethodAudit({
  hypothesis,
  metadata,
  profile,
  runMode,
  pValueThreshold,
  scoreThreshold,
  scanReverseComplement,
  completedAt,
}: {
  hypothesis: string
  metadata: MotifMetadata | null
  profile: FastaProfile | null
  runMode: RunMode
  pValueThreshold: string
  scoreThreshold: string
  scanReverseComplement: boolean
  completedAt: string
}) {
  const background = metadata?.background_frequencies

  return (
    <div className="grid gap-4 p-4 xl:grid-cols-[minmax(0,1fr)_360px]">
      <div className="space-y-4">
        <div className="rounded-[8px] border border-stone-200 bg-stone-50 p-4">
          <p className="text-xs font-semibold uppercase text-emerald-700">Scientific method</p>
          <dl className="mt-3 grid gap-3 sm:grid-cols-2">
            <AuditItem label="Hypothesis" value={hypothesis || "Not recorded"} />
            <AuditItem label="Observation unit" value={profile?.fileName ?? metadata?.filename ?? "No FASTA selected"} />
            <AuditItem label="Motif model" value={metadata?.motif_pattern ?? "No motif tested"} />
            <AuditItem
              label="Decision rule"
              value={
                runMode === "pvalue"
                  ? `P-value <= ${pValueThreshold}`
                  : `Integer score >= ${scoreThreshold || "not set"}`
              }
            />
            <AuditItem label="Reverse complement" value={scanReverseComplement ? "Included" : "Not included"} />
            <AuditItem label="Completed" value={completedAt ? new Date(completedAt).toLocaleString() : "No completed run"} />
          </dl>
        </div>

        <div className="rounded-[8px] border border-stone-200 bg-white p-4">
          <p className="text-sm font-semibold text-stone-950">Article-aligned interpretation</p>
          <div className="mt-3 grid gap-3 md:grid-cols-3">
            <Principle
              title="Background-aware"
              body="The scan estimates A/C/G/T frequencies from the uploaded FASTA, then evaluates each score against that null model."
            />
            <Principle
              title="Threshold explicit"
              body="A p-value threshold makes the decision rule auditable and avoids treating raw PWM scores as standalone evidence."
            />
            <Principle
              title="Approximation visible"
              body="The backend uses an integer log-odds PWM, so the UI keeps the chosen cutoff mode and score scale in the audit trail."
            />
          </div>
        </div>
      </div>

      <div className="rounded-[8px] border border-stone-200 bg-white p-4">
        <p className="text-sm font-semibold text-stone-950">Background frequencies</p>
        {background ? (
          <div className="mt-4 space-y-3">
            {BASES.map((base) => (
              <div className="grid grid-cols-[28px_minmax(0,1fr)_52px] items-center gap-3" key={base}>
                <span className={`flex size-7 items-center justify-center rounded-[6px] border text-xs font-bold ${BASE_SWATCHES[base]}`}>
                  {base}
                </span>
                <div className="h-2 rounded-full bg-stone-100">
                  <div
                    className="h-2 rounded-full"
                    style={{
                      backgroundColor: BASE_HEX[base],
                      width: `${Math.max(2, (background[base] ?? 0) * 100)}%`,
                    }}
                  />
                </div>
                <span className="text-right font-mono text-xs text-stone-700">
                  {formatPercent(background[base] ?? 0)}
                </span>
              </div>
            ))}
          </div>
        ) : (
          <p className="mt-4 text-sm text-stone-500">No background model has been computed yet.</p>
        )}

        {metadata ? (
          <dl className="mt-5 space-y-2 border-t border-stone-200 pt-4 text-sm">
            <AuditLine label="Motif length" value={String(metadata.motif_length)} />
            <AuditLine label="Cutoff mode" value={metadata.cutoff_mode} />
            <AuditLine label="Score threshold" value={metadata.score_threshold === null ? "Derived from p-value" : String(metadata.score_threshold)} />
          </dl>
        ) : null}
      </div>
    </div>
  )
}

function AuditItem({ label, value }: { label: string; value: string }) {
  return (
    <div>
      <dt className="text-xs font-medium uppercase text-stone-500">{label}</dt>
      <dd className="mt-1 text-sm text-stone-900">{value}</dd>
    </div>
  )
}

function AuditLine({ label, value }: { label: string; value: string }) {
  return (
    <div className="flex items-center justify-between gap-3">
      <dt className="text-stone-500">{label}</dt>
      <dd className="font-mono text-stone-900">{value}</dd>
    </div>
  )
}

function Principle({ title, body }: { title: string; body: string }) {
  return (
    <div className="rounded-[8px] border border-stone-200 bg-stone-50 p-3">
      <h3 className="text-sm font-semibold text-stone-950">{title}</h3>
      <p className="mt-2 text-sm leading-6 text-stone-600">{body}</p>
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
  const bestHit = [...matches].sort((a, b) => sortMatches(a, b, "pvalue"))[0] ?? null
  const scores = matches.map((match) => match.score)

  return {
    totalMatches: matches.length,
    uniqueGenes: new Set(matches.map((match) => match.gene)).size,
    forwardSource: matches.filter((match) => match.source_strand === "forward").length,
    reverseSource: matches.filter((match) => match.source_strand === "reverse").length,
    bestHit,
    maxScore: scores.length ? Math.max(...scores) : null,
    medianPValue: pValues.length ? pValues[Math.floor(pValues.length / 2)] : null,
  }
}

function summarizeGenes(matches: MotifMatch[]): GeneSummaryRow[] {
  const grouped = new Map<string, MotifMatch[]>()

  for (const match of matches) {
    grouped.set(match.gene, [...(grouped.get(match.gene) ?? []), match])
  }

  return [...grouped.entries()]
    .map(([gene, geneMatches]) => {
      const sorted = [...geneMatches].sort((a, b) => sortMatches(a, b, "pvalue"))
      const firstByPosition = [...geneMatches].sort((a, b) => a.source_position - b.source_position)[0]
      const sourceStrands = [...new Set(geneMatches.map((match) => match.source_strand_symbol))].join(" / ")

      return {
        gene,
        hits: geneMatches.length,
        bestPValue: sorted[0].p_value_numeric,
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

function parseScanResponse(raw: string): ApiPayload {
  const trimmed = raw.trim()
  if (!trimmed) throw new Error("The API returned an empty response.")

  try {
    return JSON.parse(trimmed) as ApiPayload
  } catch {
    const payloads = trimmed
      .split(/\r?\n/)
      .filter(Boolean)
      .map((line) => JSON.parse(line) as ApiPayload)
    const payload = payloads.find((candidate) => candidate.metadata && Array.isArray(candidate.matches))
    if (!payload) throw new Error("The API response could not be parsed.")
    return payload
  }
}

function readApiError(raw: string) {
  try {
    const parsed = JSON.parse(raw) as { detail?: string }
    return parsed.detail ?? raw
  } catch {
    return raw || "The API rejected the scan request."
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
  if (!Number.isFinite(value)) return "—"
  return value.toExponential(2)
}

function escapeCsv(value: string) {
  return `"${value.replaceAll('"', '""')}"`
}

function safeTimestamp() {
  return new Date().toISOString().replaceAll(":", "-").replace(/\.\d{3}Z$/, "Z")
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
