export class TSVector {
  constructor(private set = new Set<string>()) {}
  add = (other: string) => {
    this.set.add(other)
  }
  intersect(other: TSVector) {
    const intersection = new Set<string>()
    if (this.size < other.size) {
      this.set.forEach(el => {
        if (other.set.has(el)) intersection.add(el)
      })
    } else {
      other.set.forEach(el => {
        if (this.set.has(el)) intersection.add(el)
      })
    }
    return new TSVector(intersection)
  }
  get size() {
    return this.set.size
  }
}

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
export default function tsvector(s: string) {
  s = s.toLowerCase()
    .replace(/\s+/g, '  ')
    .replace(/^\s*/g, s.length > 1 ? ' ' : '  ')
    .replace(/\s*$/g, '  ')
  const trigrams = new TSVector()
  for (let i = 0, j = 3; j <= s.length; i++, j++) {
    const ngram = s.slice(i, j)
    if (ngram.endsWith('  ')) continue
    trigrams.add(ngram)
  }
  return trigrams
}