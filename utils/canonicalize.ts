import { canonicalize as canonicalizeRFC8785 } from 'json-canonicalize'

/**
 * Mirror JSON.stringify's treatment of non-finite numbers (NaN/Infinity -> null).
 *
 * json-canonicalize >=3 follows RFC 8785 strictly and throws on non-finite numbers,
 *  earlier versions serialized them as null. Content-addressed ids are derived from
 *  this output, so we keep the earlier behavior: results (e.g. statistics) may well
 *  contain NaN, and the ids of existing records must not change.
 */
function replaceNonFinite(value: unknown): unknown {
  if (typeof value === 'number') return Number.isFinite(value) ? value : null
  if (value === null || typeof value !== 'object') return value
  // like JSON.stringify, toJSON (e.g. on Date) takes precedence -- but json-canonicalize
  //  serializes toJSON results differently (without sorting keys), so we keep a toJSON
  //  object rather than returning the plain result which would change the output
  if (typeof (value as { toJSON?: unknown }).toJSON === 'function') {
    const json = replaceNonFinite((value as { toJSON: () => unknown }).toJSON())
    return { toJSON: () => json }
  }
  if (Array.isArray(value)) return value.map(replaceNonFinite)
  const result: Record<string, unknown> = {}
  for (const key of Object.keys(value)) {
    result[key] = replaceNonFinite((value as Record<string, unknown>)[key])
  }
  return result
}

/**
 * RFC 8785 JSON canonicalization, with NaN/Infinity serialized as null
 */
export function canonicalize(data: unknown): string {
  return canonicalizeRFC8785(replaceNonFinite(data))
}
