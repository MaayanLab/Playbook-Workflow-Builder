/**
 * tri-gram searching works as follows:
 * Given: "my text"
 *  we produce tokens like so: "  my  text " => {"  m", " my", "my ", "  t", " te", ...}
 * Then we do the same to the search and compute a set overlap.
 * Two spaces are used in front of a word boundary to capture single letter prefix matching,
 *   but two spaces are not used for the postfix of a word boundary.
 * This variant of trigram searching seems to produce the most intuitive results when using
 *  1 or 2 characters for the search.
 */
export function tsvector(s: string) {
  s = s.toLowerCase()
    .replace(/\s+/g, '  ')
    .replace(/^\s*/g, s.length > 1 ? ' ' : '  ')
    .replace(/\s*$/g, '  ')
  const trigrams = new Set<string>()
  for (let i = 0, j = 3; j <= s.length; i++, j++) {
    const ngram = s.slice(i, j)
    if (ngram.endsWith('  ')) continue
    trigrams.add(ngram)
  }
  return trigrams
}

export function tsvector_intersect(a: Set<string>, b: Set<string>) {
  const intersection = new Set<string>()
  // choose the smaller index to loop through
  if (a.size < b.size) {
    // collect all intersecting trigrams
    a.forEach(el => {
      if (b.has(el)) intersection.add(el)
    })
  } else {
    // collect all intersecting trigrams
    b.forEach(el => {
      if (a.has(el)) intersection.add(el)
    })
  }
  return intersection
}
