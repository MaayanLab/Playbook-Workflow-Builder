function matchAll(str: string, re: RegExp) {
  const matches = []
  let match
  while ((match = re.exec(str)) != null) {
    matches.push(match)
  }
  return matches
}

export default function Linkify({ children }: { children: string }) {
  if (!children) return children
  const { i, el } = matchAll(children, /(doi:|https:\/\/)([^\s]+)/g)
    .reduce(({ i, el }, m) => m.index ? ({
      i: m.index + m[0].length,
      el: [
        ...el,
        children.slice(i, m.index),
        m[1] === 'doi:' ? <a href={`https://doi.org/${m[2]}`}>{m[0]}</a>
        : m[1] === 'http://' ? <a href={m[0]}>{m[0]}</a>
        : m[1] === 'https://' ? <a href={m[0]}>{m[0]}</a>
        : m[0],
      ]
    }) : { i, el }, { i: 0, el: [] as React.ReactNode[] })
  el.push(children.slice(i))
  return el
}
