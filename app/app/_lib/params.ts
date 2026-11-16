/**
 * Route params as received by pages are still URL-encoded (e.g. `Term%5BGene%5D`),
 *  unlike route handler params and the pages router's `router.query`. Specs contain
 *  characters like `[`, `]` and spaces, so pages must decode them before use.
 */
export function decodeParam(value: string) {
  try {
    return decodeURIComponent(value)
  } catch {
    // not a valid encoding (e.g. a literal `%`), use it as is
    return value
  }
}
