import type { AnalysisDetail, Hit } from "@/types/api"

const API_BASE = import.meta.env.VITE_API_BASE_URL ?? "http://127.0.0.1:8000/api"

async function parseResponse<T>(response: Response): Promise<T> {
  if (!response.ok) {
    const body = await response.json().catch(() => null)
    const detail = body?.detail ?? `HTTP ${response.status}`
    throw new Error(detail)
  }
  return response.json() as Promise<T>
}

export async function createAnalysis(
  file: File,
  motif: string,
  pValueCutoff?: number,
): Promise<AnalysisDetail> {
  const payload = new FormData()
  payload.append("file", file)
  payload.append("motif", motif)
  if (typeof pValueCutoff === "number" && Number.isFinite(pValueCutoff)) {
    payload.append("p_value_cutoff", String(pValueCutoff))
  }

  const response = await fetch(`${API_BASE}/analyses`, {
    method: "POST",
    body: payload,
  })
  return parseResponse<AnalysisDetail>(response)
}

export async function fetchResults(analysisId: string, limit = 100, offset = 0): Promise<Hit[]> {
  const params = new URLSearchParams({ limit: String(limit), offset: String(offset) })
  const response = await fetch(`${API_BASE}/analyses/${analysisId}/results?${params}`)
  return parseResponse<Hit[]>(response)
}
