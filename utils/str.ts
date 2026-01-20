export function truncate(s: string, length: number) {
  if (s.length <= length) return s
  else return `${s.slice(0, length-3)}...`
}
