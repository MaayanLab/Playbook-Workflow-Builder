class TSVector {
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
export default function tsvector(s: string) {
  const s_split = ` ${s.toLowerCase().replace(/^\s+/g, '').replace(/\s+$/g, '').replace(/\s+/g, ' ')} `.split('')
  const trigrams = new TSVector()
  for (let i = 0, j = 3; j < s.length; i++, j++) {
    trigrams.add(s_split.slice(i, j).join(''))
  }
  return trigrams
}