import { useState, type ChangeEvent, type DragEvent } from "react"
import { zodResolver } from "@hookform/resolvers/zod"
import { useMutation } from "@tanstack/react-query"
import { useForm } from "react-hook-form"
import { z } from "zod"
import {
  CircleHelp,
  Dna,
  Loader2,
  Play,
  Upload,
} from "lucide-react"

import { createAnalysis } from "@/api/client"
import { Button } from "@/components/ui/button"
import {
  Card,
  CardContent,
  CardDescription,
  CardHeader,
  CardTitle,
} from "@/components/ui/card"
import { FieldError } from "@/components/ui/field"
import { Input } from "@/components/ui/input"
import { Label } from "@/components/ui/label"
import { Skeleton } from "@/components/ui/skeleton"
import {
  Table,
  TableBody,
  TableCell,
  TableHead,
  TableHeader,
  TableRow,
} from "@/components/ui/table"
import {
  Tooltip,
  TooltipContent,
  TooltipProvider,
  TooltipTrigger,
} from "@/components/ui/tooltip"
import type { AnalysisDetail, Hit } from "@/types/api"
import { formatFimoPValue, formatNumber, formatScientific } from "@/utils/format"

const iupacPattern = /^[ACGTRYSWKMBDHVN]+$/i

function isSupportedFastaFile(name: string) {
  const lowerName = name.toLowerCase()
  return (
    lowerName.endsWith(".fasta")
    || lowerName.endsWith(".fa")
    || lowerName.endsWith(".fna")
    || lowerName.endsWith(".fsa")
    || lowerName.endsWith(".fas")
    || lowerName.endsWith(".txt")
  )
}

function formatFileSize(sizeInBytes: number) {
  if (!Number.isFinite(sizeInBytes) || sizeInBytes < 0) {
    return "Unknown size"
  }
  if (sizeInBytes < 1024) {
    return `${sizeInBytes} B`
  }
  if (sizeInBytes < 1024 ** 2) {
    return `${(sizeInBytes / 1024).toFixed(1)} KB`
  }
  return `${(sizeInBytes / 1024 ** 2).toFixed(2)} MB`
}

const analysisFormSchema = z.object({
  file: z
    .string()
    .trim()
    .min(1, "Select a FASTA file before running.")
    .refine(
      (value) => isSupportedFastaFile(value),
      "Use a FASTA file (.fasta, .fa, .fna, .fsa, .fas) or .txt.",
    ),
  motif: z
    .string()
    .trim()
    .min(1, "Enter an IUPAC motif before running.")
    .regex(iupacPattern, "Use only IUPAC DNA letters (ACGTRYSWKMBDHVN)."),
  pValueCutoff: z
    .string()
    .trim()
    .min(1, "Enter a p-value cutoff.")
    .refine((value) => {
      const parsed = Number(value)
      return Number.isFinite(parsed) && parsed > 0 && parsed <= 1
    }, "P-value cutoff must be a number greater than 0 and up to 1."),
})

type AnalysisFormValues = z.infer<typeof analysisFormSchema>

const fieldHelpText = {
  fasta: "FASTA with promoter/upstream sequences. One sequence per record.",
  motif: "IUPAC consensus used to derive the PWM for scanning.",
  pValueCutoff: "Maximum per-occurrence p-value accepted as a motif hit.",
} as const

const summaryHelpText = {
  sequences: "Number of sequences scanned in the uploaded FASTA.",
  observedHits: "Total motif occurrences passing the p-value cutoff.",
  enrichmentPValue: "One-sided p-value for enrichment (observed > expected).",
  enrichmentQValue: "FDR-adjusted enrichment p-value across motif tests.",
  depletionPValue: "One-sided p-value for depletion (observed < expected).",
  twoSidedPValue: "Two-sided significance combining enrichment and depletion tails.",
  foldChange: "Observed/Expected hit ratio. <1 suggests depletion, >1 enrichment.",
} as const

const tableHelpText = {
  gene: "Gene identifier from the sequence record.",
  strand: "Match strand relative to the scanned sequence (+ or -).",
  position: "Start and end coordinates of the match in the scanned sequence.",
  globalLocus: "Genome coordinates derived from hit local position and the sequence genomic interval in description.",
  sequence: "Matched sequence for the reported motif occurrence.",
  pValue: "Per-occurrence p-value under the background model.",
  qValue: "FDR-adjusted per-occurrence p-value.",
} as const

type GlobalLocus = {
  contig: string | null
  regionStart: number
  regionEnd: number
  globalStart: number
  globalEnd: number
  barLeftPercent: number
  barWidthPercent: number
}

function clamp(value: number, minValue: number, maxValue: number) {
  return Math.min(Math.max(value, minValue), maxValue)
}

function parseGlobalLocus(hit: Hit): GlobalLocus | null {
  const rangeMatch = hit.description.match(/\[(\d+)-(\d+)\]/)
  if (!rangeMatch) {
    return null
  }

  const parsedStart = Number.parseInt(rangeMatch[1], 10)
  const parsedEnd = Number.parseInt(rangeMatch[2], 10)
  if (!Number.isFinite(parsedStart) || !Number.isFinite(parsedEnd)) {
    return null
  }

  const regionStart = Math.min(parsedStart, parsedEnd)
  const regionEnd = Math.max(parsedStart, parsedEnd)
  const regionLength = Math.max(regionEnd - regionStart + 1, 1)
  const hitLocalStart = Math.min(hit.start, hit.end)
  const hitLocalEnd = Math.max(hit.start, hit.end)
  const isReverseComplementRegion = /reverse complement/i.test(hit.description)

  const mappedStart = isReverseComplementRegion
    ? regionEnd - hitLocalEnd + 1
    : regionStart + hitLocalStart - 1
  const mappedEnd = isReverseComplementRegion
    ? regionEnd - hitLocalStart + 1
    : regionStart + hitLocalEnd - 1

  const globalStart = Math.min(mappedStart, mappedEnd)
  const globalEnd = Math.max(mappedStart, mappedEnd)
  const highlightLength = Math.max(globalEnd - globalStart + 1, 1)

  const rawLeftPercent = ((globalStart - regionStart) / regionLength) * 100
  const leftPercent = clamp(rawLeftPercent, 0, 100)
  const rawWidthPercent = (highlightLength / regionLength) * 100
  const maxAvailableWidth = Math.max(100 - leftPercent, 0)
  const widthPercent = clamp(rawWidthPercent, Math.min(0.8, maxAvailableWidth), maxAvailableWidth)

  const contigMatch = hit.description.match(/\b(supercont[^\s\]]+|NW_[\w.]+|chr[^\s\]]+|scaffold[^\s\]]+|contig[^\s\]]+)\b/i)

  return {
    contig: contigMatch?.[1] ?? null,
    regionStart,
    regionEnd,
    globalStart,
    globalEnd,
    barLeftPercent: leftPercent,
    barWidthPercent: widthPercent,
  }
}

export function MotifScannerPage() {
  const [isDragActive, setIsDragActive] = useState(false)
  const [selectedFile, setSelectedFileState] = useState<File | null>(null)
  const [analysis, setAnalysis] = useState<AnalysisDetail | null>(null)

  const {
    register,
    setValue,
    handleSubmit,
    clearErrors,
    setError,
    watch,
    formState: { errors },
  } = useForm<AnalysisFormValues>({
    resolver: zodResolver(analysisFormSchema),
    defaultValues: {
      file: "",
      motif: "",
      pValueCutoff: "1e-4",
    },
  })

  const mutation = useMutation({
    mutationFn: async (values: AnalysisFormValues) => {
      const targetFile = selectedFile
      if (!targetFile) {
        throw new Error("Select a FASTA file before running.")
      }
      const targetMotif = values.motif.trim().toUpperCase()
      const pValueCutoff = Number(values.pValueCutoff)
      return createAnalysis(targetFile, targetMotif, pValueCutoff)
    },
    onMutate: () => {
      setAnalysis(null)
    },
    onSuccess: (result) => {
      setAnalysis(result)
    },
    onError: (requestError) => {
      setError("root.serverError", {
        type: "server",
        message: requestError instanceof Error ? requestError.message : "Failed to create analysis.",
      })
    },
  })

  const motifField = register("motif", {
    onChange: () => clearErrors("root.serverError"),
  })
  const pValueCutoffField = register("pValueCutoff", {
    onChange: () => clearErrors("root.serverError"),
  })

  function setSelectedFile(nextFile: File | null) {
    const safeFile = nextFile instanceof File ? nextFile : null
    setValue(
      "file",
      safeFile?.name ?? "",
      {
        shouldDirty: true,
        shouldTouch: true,
        shouldValidate: true,
      },
    )
    setSelectedFileState(safeFile)
    clearErrors("root.serverError")
  }

  function handleFileInputChange(event: ChangeEvent<HTMLInputElement>) {
    setSelectedFile(event.target.files?.[0] ?? null)
  }

  function handleFileDrop(event: DragEvent<HTMLLabelElement>) {
    event.preventDefault()
    setIsDragActive(false)
    setSelectedFile(event.dataTransfer.files?.[0] ?? null)
  }

  function handleFileDragOver(event: DragEvent<HTMLLabelElement>) {
    event.preventDefault()
    event.dataTransfer.dropEffect = "copy"
    setIsDragActive(true)
  }

  function handleFileDragLeave(event: DragEvent<HTMLLabelElement>) {
    event.preventDefault()
    setIsDragActive(false)
  }

  function normalizeCutoffInput() {
    const parsed = Number(watch("pValueCutoff"))
    if (Number.isFinite(parsed) && parsed > 0 && parsed <= 1) {
      setValue("pValueCutoff", parsed.toExponential().replace("e+", "e"), {
        shouldValidate: true,
      })
    }
  }

  function onSubmit(values: AnalysisFormValues) {
    clearErrors("root.serverError")
    mutation.mutate(values)
  }

  return (
    <TooltipProvider delayDuration={150}>
      <main className="min-h-screen bg-background">
        <div className="mx-auto flex w-full max-w-7xl flex-col gap-6 px-5 py-6 lg:px-8">
          <header className="flex flex-col gap-4 border-b pb-5 lg:flex-row lg:items-end lg:justify-between">
            <div className="max-w-3xl">
              <div className="mb-2 flex items-center gap-2 text-sm font-medium text-muted-foreground">
                <Dna className="h-4 w-4 text-primary" />
                Motif Scan
              </div>
              <h1 className="text-2xl font-semibold tracking-normal text-foreground">
                Statistical Motif Search in Promoter Regions
              </h1>
              <p className="mt-2 text-sm leading-6 text-muted-foreground">
                Upload a FASTA file, enter an IUPAC motif, and rank hits by p-value,
                q-value, and enrichment.
              </p>
            </div>
          </header>

          <section className="grid gap-6 lg:grid-cols-[400px_1fr] lg:items-start">
            <Card className="h-min self-start">
              <CardHeader>
                <CardTitle>Parameters</CardTitle>
                <CardDescription>FASTA, IUPAC motif, and p-value cutoff.</CardDescription>
              </CardHeader>
              <CardContent className="space-y-5">
                <form className="space-y-5" onSubmit={handleSubmit(onSubmit)} noValidate>
                  <div className="space-y-2">
                    <input type="hidden" {...register("file")} />
                    <div className="flex items-center gap-1.5">
                      <Label htmlFor="fasta">FASTA File</Label>
                      <InlineHelp text={fieldHelpText.fasta} />
                    </div>
                    <label
                      htmlFor="fasta"
                      onDrop={handleFileDrop}
                      onDragEnter={handleFileDragOver}
                      onDragOver={handleFileDragOver}
                      onDragLeave={handleFileDragLeave}
                      className={`flex w-full min-w-0 cursor-pointer flex-col items-center justify-center gap-2 rounded-md border border-dashed p-5 text-center transition-colors ${
                        errors.file
                          ? "border-destructive bg-destructive/5"
                          : isDragActive
                            ? "border-primary bg-primary/10"
                            : "border-border bg-muted/30 hover:border-primary/60 hover:bg-muted/60"
                      }`}
                    >
                      <Input
                        id="fasta"
                        type="file"
                        accept=".fasta,.fa,.fna,.fsa,.fas,.txt"
                        className="sr-only"
                        onChange={handleFileInputChange}
                      />
                      <Upload className="h-5 w-5 text-primary" />
                      <div className="w-full min-w-0 space-y-1">
                        <p
                          className="truncate text-sm font-medium text-foreground"
                          title={selectedFile?.name ?? undefined}
                        >
                          {selectedFile ? selectedFile.name : "Drop FASTA here or click to upload"}
                        </p>
                        <p className="text-xs text-muted-foreground">
                          {selectedFile
                            ? `${formatFileSize(selectedFile.size)} - ${selectedFile.type || "local file"}`
                            : "Accepted: .fasta, .fa, .fna, .fsa, .fas, .txt"}
                        </p>
                      </div>
                    </label>
                    <FieldError errors={[errors.file]} />
                  </div>

                  <div className="space-y-2">
                    <div className="flex items-center gap-1.5">
                      <Label htmlFor="motif">IUPAC Motif</Label>
                      <InlineHelp text={fieldHelpText.motif} />
                    </div>
                    <Input
                      id="motif"
                      spellCheck={false}
                      placeholder="IUPAC Motif"
                      className="font-mono uppercase"
                      {...motifField}
                    />
                    <FieldError errors={[errors.motif]} />
                  </div>

                  <div className="space-y-2">
                    <div className="flex items-center gap-1.5">
                      <Label htmlFor="pvalue-cutoff">P-value cutoff</Label>
                      <InlineHelp text={fieldHelpText.pValueCutoff} />
                    </div>
                    <Input
                      id="pvalue-cutoff"
                      type="text"
                      inputMode="decimal"
                      className="font-mono"
                      placeholder="1e-4"
                      {...pValueCutoffField}
                      onBlur={(event) => {
                        pValueCutoffField.onBlur(event)
                        normalizeCutoffInput()
                      }}
                    />
                    <FieldError errors={[errors.pValueCutoff]} />
                  </div>

                  <div>
                    <Button
                      className="w-full justify-center font-semibold"
                      size="lg"
                      type="submit"
                      disabled={mutation.isPending}
                    >
                      {mutation.isPending ? (
                        <Loader2 className="mr-2 h-4 w-4 animate-spin" />
                      ) : (
                        <Play className="mr-2 h-4 w-4" />
                      )}
                      {mutation.isPending ? "Scanning..." : "Run Scan"}
                    </Button>
                    <FieldError errors={[errors.root?.serverError]} className="mt-2" />
                  </div>
                </form>
              </CardContent>
            </Card>

            <div className="space-y-6">
              <SummaryCards analysis={analysis} loading={mutation.isPending} />
              <HitsTable analysis={analysis} loading={mutation.isPending} />
            </div>
          </section>
        </div>
      </main>
    </TooltipProvider>
  )
}

function SummaryCards({ analysis, loading }: { analysis: AnalysisDetail | null; loading: boolean }) {
  const summary = analysis?.summary

  const metrics = [
    {
      label: "Sequences",
      value: summary ? formatNumber(summary.sequence_count, 0) : "-",
      help: summaryHelpText.sequences,
    },
    {
      label: "Observed Hits",
      value: loading ? "..." : summary ? formatNumber(summary.observed_hits, 0) : "-",
      help: summaryHelpText.observedHits,
    },
    {
      label: "Enrichment p-value",
      value: summary ? formatScientific(summary.enrichment_p_value) : "-",
      help: summaryHelpText.enrichmentPValue,
    },
    {
      label: "Enrichment q-value",
      value: summary ? formatScientific(summary.enrichment_q_value) : "-",
      help: summaryHelpText.enrichmentQValue,
    },
    {
      label: "Depletion p-value",
      value: summary ? formatScientific(summary.depletion_p_value) : "-",
      help: summaryHelpText.depletionPValue,
    },
    {
      label: "Two-sided p-value",
      value: summary ? formatScientific(summary.two_sided_p_value) : "-",
      help: summaryHelpText.twoSidedPValue,
    },
    {
      label: "Fold-change (Obs/Exp)",
      value: summary ? formatNumber(summary.fold_change, 3) : "-",
      help: summaryHelpText.foldChange,
    },
  ]

  return (
    <Card>
      <CardContent className="p-3 sm:p-4">
        {loading ? (
          <div className="space-y-2">
            <Skeleton className="h-8 w-full" />
            <Skeleton className="h-8 w-10/12" />
          </div>
        ) : (
          <div className="grid gap-2 sm:grid-cols-2 xl:grid-cols-3">
            {metrics.map((metric) => (
              <div
                key={metric.label}
                className="flex items-center justify-between gap-3 rounded-sm border bg-muted/20 px-3 py-2"
              >
                <div className="flex min-w-0 items-center gap-1.5">
                  <p className="truncate text-xs text-muted-foreground">{metric.label}</p>
                  <InlineHelp text={metric.help} />
                </div>
                <p className="truncate text-sm font-semibold text-foreground" title={metric.value}>
                  {metric.value}
                </p>
              </div>
            ))}
          </div>
        )}
      </CardContent>
    </Card>
  )
}

function HitsTable({ analysis, loading }: { analysis: AnalysisDetail | null; loading: boolean }) {
  if (loading) {
    return (
      <Card>
        <CardHeader>
          <CardTitle>Prioritized PWM Hits</CardTitle>
          <CardDescription>Loading new scan results...</CardDescription>
        </CardHeader>
        <CardContent>
          <div className="space-y-3">
            <Skeleton className="h-6 w-52" />
            <Skeleton className="h-[320px] w-full" />
          </div>
        </CardContent>
      </Card>
    )
  }

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
    const pValueDelta = left.window_p_value - right.window_p_value
    if (pValueDelta !== 0) {
      return pValueDelta
    }
    return left.q_value - right.q_value
  })

  return (
    <Card>
      <CardHeader>
        <CardTitle>Prioritized PWM Hits</CardTitle>
        <CardDescription>Sequences sorted by lowest p-value.</CardDescription>
      </CardHeader>
      <CardContent>
        <div className="h-[430px] overflow-auto rounded-md border">
          <Table>
            <TableHeader>
              <TableRow>
                <TableHead className="sticky top-0 z-10 bg-background">
                  <HeaderWithHelp label="gene" help={tableHelpText.gene} />
                </TableHead>
                <TableHead className="sticky top-0 z-10 bg-background">
                  <HeaderWithHelp label="strand" help={tableHelpText.strand} />
                </TableHead>
                <TableHead className="sticky top-0 z-10 bg-background">
                  <HeaderWithHelp label="position" help={tableHelpText.position} />
                </TableHead>
                <TableHead className="sticky top-0 z-10 bg-background">
                  <HeaderWithHelp label="global locus" help={tableHelpText.globalLocus} />
                </TableHead>
                <TableHead className="sticky top-0 z-10 bg-background">
                  <HeaderWithHelp label="sequence" help={tableHelpText.sequence} />
                </TableHead>
                <TableHead className="sticky top-0 z-10 bg-background">
                  <HeaderWithHelp label="p-value" help={tableHelpText.pValue} />
                </TableHead>
                <TableHead className="sticky top-0 z-10 bg-background">
                  <HeaderWithHelp label="q-value" help={tableHelpText.qValue} />
                </TableHead>
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
                  <TableCell className="min-w-[260px]">
                    <GlobalLocusBar hit={hit} />
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

function HeaderWithHelp({ label, help }: { label: string; help: string }) {
  return (
    <div className="inline-flex items-center gap-1.5">
      <span>{label}</span>
      <InlineHelp text={help} />
    </div>
  )
}

function InlineHelp({ text }: { text: string }) {
  return (
    <Tooltip>
      <TooltipTrigger asChild>
        <button
          type="button"
          className="inline-flex h-4 w-4 items-center justify-center text-muted-foreground/80 transition-colors hover:text-foreground"
          aria-label="Show info"
        >
          <CircleHelp className="h-3.5 w-3.5" />
        </button>
      </TooltipTrigger>
      <TooltipContent side="top" sideOffset={6} className="max-w-[280px] leading-relaxed">
        {text}
      </TooltipContent>
    </Tooltip>
  )
}

function GlobalLocusBar({ hit }: { hit: Hit }) {
  const locus = parseGlobalLocus(hit)
  if (!locus) {
    return <span className="text-xs text-muted-foreground">N/A</span>
  }

  const globalRangeLabel = `${locus.globalStart}-${locus.globalEnd}`
  const regionLabel = `${locus.regionStart}-${locus.regionEnd}`

  return (
    <div className="space-y-1.5">
      <p className="truncate text-xs text-foreground" title={`${locus.contig ?? "contig?"}:${globalRangeLabel}`}>
        {(locus.contig ?? "contig?")}:{globalRangeLabel}
      </p>
      <div className="space-y-1">
        <div className="relative h-2 w-full overflow-hidden rounded-sm bg-muted">
          <div
            className="absolute top-0 h-full rounded-sm bg-blue-500"
            style={{
              left: `${locus.barLeftPercent}%`,
              width: `${locus.barWidthPercent}%`,
            }}
          />
        </div>
        <p className="truncate text-[11px] text-muted-foreground" title={`Region: ${regionLabel}`}>
          Region: {regionLabel}
        </p>
      </div>
    </div>
  )
}
