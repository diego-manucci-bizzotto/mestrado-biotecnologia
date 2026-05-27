export function formatNumber(value: number, digits = 3): string {
  return new Intl.NumberFormat("pt-BR", {
    maximumFractionDigits: digits,
  }).format(value)
}

export function formatScientific(value: number): string {
  if (value < 0.01 || value >= 10_000) {
    return value.toExponential(2)
  }
  return formatNumber(value, 6)
}

export function formatFimoPValue(value: number): string {
  const scientific = value.toExponential(2)
  const [mantissa, exponentRaw] = scientific.split("e")
  const sign = exponentRaw.startsWith("-") ? "-" : "+"
  const exponentDigits = exponentRaw.replace(/[+-]/g, "").padStart(2, "0")
  return `${mantissa}e${sign}${exponentDigits}`
}
