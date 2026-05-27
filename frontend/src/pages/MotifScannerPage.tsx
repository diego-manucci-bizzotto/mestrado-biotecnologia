import { useState } from "react"
import {
  Activity,
  AlertCircle,
  BadgeCheck,
  Database,
  Dna,
  Loader2,
  Upload,
} from "lucide-react"

import { createAnalysis } from "@/api/client"
import { Alert, AlertDescription, AlertTitle } from "@/components/ui/alert"
import { Button } from "@/components/ui/button"
import {
  Card,
  CardContent,
  CardDescription,
  CardHeader,
  CardTitle,
} from "@/components/ui/card"
import { Input } from "@/components/ui/input"
import { Label } from "@/components/ui/label"
import {
  Table,
  TableBody,
  TableCell,
  TableHead,
  TableHeader,
  TableRow,
} from "@/components/ui/table"
import type { AnalysisDetail } from "@/types/api"
import { formatFimoPValue, formatNumber, formatScientific } from "@/utils/format"

export function MotifScannerPage() {
  const [motifInput, setMotifInput] = useState("")
  const [pValueCutoffInput, setPValueCutoffInput] = useState("1e-4")
  const [file, setFile] = useState<File | null>(null)
  const [analysis, setAnalysis] = useState<AnalysisDetail | null>(null)
  const [loading, setLoading] = useState(false)
  const [error, setError] = useState<string | null>(null)

  async function runAnalysis(
    targetFile = file,
    targetMotif = motifInput,
  ) {
    if (!targetFile) {
      setError("Select a FASTA file before running.")
      return
    }
    if (!targetMotif.trim()) {
      setError("Enter an IUPAC motif before running.")
      return
    }
    const pValueCutoff = Number(pValueCutoffInput)
    if (!Number.isFinite(pValueCutoff) || pValueCutoff <= 0 || pValueCutoff > 1) {
      setError("P-value cutoff must be a number greater than 0 and up to 1.")
      return
    }

    setLoading(true)
    setError(null)
    try {
      const result = await createAnalysis(targetFile, targetMotif, pValueCutoff)
      setAnalysis(result)
    } catch (requestError) {
      setError(requestError instanceof Error ? requestError.message : "Failed to create analysis.")
    } finally {
      setLoading(false)
    }
  }

  function normalizeCutoffInput() {
    const parsed = Number(pValueCutoffInput)
    if (Number.isFinite(parsed) && parsed > 0 && parsed <= 1) {
      setPValueCutoffInput(parsed.toExponential().replace("e+", "e"))
    }
  }

  return (
    <main className="min-h-screen bg-background">
      <div className="mx-auto flex w-full max-w-7xl flex-col gap-6 px-5 py-6 lg:px-8">
        <header className="flex flex-col gap-4 border-b pb-5 lg:flex-row lg:items-end lg:justify-between">
          <div className="max-w-3xl">
            <div className="mb-2 flex items-center gap-2 text-sm font-medium text-muted-foreground">
              <Dna className="h-4 w-4 text-primary" />
              Biotechnology MSc
            </div>
            <h1 className="text-2xl font-semibold tracking-normal text-foreground">
              Statistical Motif Search in Promoter Regions
            </h1>
            <p className="mt-2 text-sm leading-6 text-muted-foreground">
              Enter an IUPAC consensus, generate a derived PWM, and prioritize matches
              by p-value, q-value, and enrichment.
            </p>
          </div>
        </header>

        {error ? (
          <Alert variant="destructive">
            <AlertCircle className="h-4 w-4" />
            <AlertTitle>Analysis Failed</AlertTitle>
            <AlertDescription>{error}</AlertDescription>
          </Alert>
        ) : null}

        <section className="grid gap-6 lg:grid-cols-[360px_1fr] lg:items-start">
          <Card className="h-min self-start">
            <CardHeader>
              <CardTitle>Parameters</CardTitle>
              <CardDescription>FASTA file and IUPAC consensus used to build the PWM.</CardDescription>
            </CardHeader>
            <CardContent className="space-y-5">
              <div className="space-y-2">
                <Label htmlFor="fasta">FASTA File</Label>
                <Input
                  id="fasta"
                  type="file"
                  accept=".fasta,.fa,.txt"
                  onChange={(event) => setFile(event.target.files?.[0] ?? null)}
                />
              </div>

              <div className="space-y-2">
                <Label htmlFor="motif">IUPAC Motif</Label>
                <Input
                  id="motif"
                  value={motifInput}
                  placeholder="e.g., GCCARG, GNGGCKCA, WCGCGWNM"
                  spellCheck={false}
                  className="font-mono uppercase"
                  onChange={(event) => setMotifInput(event.target.value.toUpperCase())}
                />
              </div>

              <div className="space-y-2">
                <Label htmlFor="pvalue-cutoff">P-value cutoff</Label>
                <Input
                  id="pvalue-cutoff"
                  type="text"
                  inputMode="decimal"
                  className="font-mono"
                  placeholder="e.g., 1e-4"
                  value={pValueCutoffInput}
                  onChange={(event) => setPValueCutoffInput(event.target.value)}
                  onBlur={normalizeCutoffInput}
                />
              </div>

              <div>
                <Button className="w-full" onClick={() => void runAnalysis()} disabled={loading}>
                  {loading ? (
                    <Loader2 className="mr-2 h-4 w-4 animate-spin" />
                  ) : (
                    <Upload className="mr-2 h-4 w-4" />
                  )}
                  Run
                </Button>
              </div>
            </CardContent>
          </Card>

          <div className="space-y-6">
            <SummaryCards analysis={analysis} loading={loading} />

            <HitsTable analysis={analysis} />
          </div>
        </section>
      </div>
    </main>
  )
}

function SummaryCards({ analysis, loading }: { analysis: AnalysisDetail | null; loading: boolean }) {
  const summary = analysis?.summary
  return (
    <div className="grid gap-3 sm:grid-cols-2 xl:grid-cols-4">
      <MetricCard
        icon={<Database className="h-4 w-4" />}
        label="Sequences"
        value={summary ? formatNumber(summary.sequence_count, 0) : "-"}
      />
      <MetricCard
        icon={<Activity className="h-4 w-4" />}
        label="Observed Hits"
        value={loading ? "..." : summary ? formatNumber(summary.observed_hits, 0) : "-"}
      />
      <MetricCard
        icon={<BadgeCheck className="h-4 w-4" />}
        label="Enrichment p-value"
        value={summary ? formatScientific(summary.enrichment_p_value) : "-"}
      />
      <MetricCard
        icon={<BadgeCheck className="h-4 w-4" />}
        label="Enrichment q-value"
        value={summary ? formatScientific(summary.enrichment_q_value) : "-"}
      />
      <MetricCard
        icon={<BadgeCheck className="h-4 w-4" />}
        label="Depletion p-value"
        value={summary ? formatScientific(summary.depletion_p_value) : "-"}
      />
      <MetricCard
        icon={<BadgeCheck className="h-4 w-4" />}
        label="Two-sided p-value"
        value={summary ? formatScientific(summary.two_sided_p_value) : "-"}
      />
      <MetricCard
        icon={<Activity className="h-4 w-4" />}
        label="Fold-change (Obs/Exp)"
        value={summary ? formatNumber(summary.fold_change, 3) : "-"}
      />
    </div>
  )
}

function MetricCard({ icon, label, value }: { icon: React.ReactNode; label: string; value: string }) {
  return (
    <Card>
      <CardContent className="flex items-center gap-3 p-4">
        <div className="rounded-md bg-primary/10 p-2 text-primary">{icon}</div>
        <div className="min-w-0">
          <p className="text-xs font-medium text-muted-foreground">{label}</p>
          <p className="truncate text-lg font-semibold">{value}</p>
        </div>
      </CardContent>
    </Card>
  )
}

function HitsTable({ analysis }: { analysis: AnalysisDetail | null }) {
  if (!analysis) {
    return (
      <Card>
        <CardContent className="p-6 text-sm text-muted-foreground">
          Run an analysis to view prioritized PWM hits.
        </CardContent>
      </Card>
    )
  }

  const sortedHits = [...analysis.hits].sort((left, right) => {
    const qValueDelta = left.q_value - right.q_value
    if (qValueDelta !== 0) {
      return qValueDelta
    }
    return left.window_p_value - right.window_p_value
  })

  return (
    <Card>
      <CardHeader>
        <CardTitle>Prioritized PWM Hits</CardTitle>
        <CardDescription>
          Sequences sorted by lowest q-value.
        </CardDescription>
      </CardHeader>
      <CardContent>
        <div className="h-[430px] overflow-y-auto rounded-md border">
          <Table>
            <TableHeader>
              <TableRow className="hover:bg-transparent">
                <TableHead className="sticky top-0 z-10 bg-background shadow-sm">Gene</TableHead>
                <TableHead className="sticky top-0 z-10 bg-background shadow-sm">Strand</TableHead>
                <TableHead className="sticky top-0 z-10 bg-background shadow-sm">Position</TableHead>
                <TableHead className="sticky top-0 z-10 bg-background shadow-sm">Sequence</TableHead>
                <TableHead className="sticky top-0 z-10 bg-background shadow-sm">P</TableHead>
                <TableHead className="sticky top-0 z-10 bg-background shadow-sm">Q</TableHead>
              </TableRow>
            </TableHeader>
            <TableBody>
              {sortedHits.slice(0, 100).map((hit) => (
                <TableRow key={hit.id}>
                  <TableCell className="font-medium">{hit.gene_id}</TableCell>
                  <TableCell>{hit.strand}</TableCell>
                  <TableCell>
                    {hit.start}-{hit.end}
                  </TableCell>
                  <TableCell>
                    <code>{hit.matched_sequence}</code>
                  </TableCell>
                  <TableCell>{formatFimoPValue(hit.window_p_value)}</TableCell>
                  <TableCell>{formatFimoPValue(hit.q_value)}</TableCell>
                </TableRow>
              ))}
            </TableBody>
          </Table>
        </div>
      </CardContent>
    </Card>
  )
}
