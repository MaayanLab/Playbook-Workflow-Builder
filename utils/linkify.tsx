import React from 'react'

function matchAll(str: string, re: RegExp) {
  const matches = []
  let match
  while ((match = re.exec(str)) != null) {
    matches.push(match)
  }
  return matches
}

export default function Linkify({ children }: { children?: string | null }) {
  if (!children) return children
  const { i, el } = matchAll(children, /(doi:|https?:\/\/)([^\s]+)/g)
    .reduce(({ i, el }, m) => m.index ? ({
      i: m.index + m[0].length,
      el: [
        ...el,
        <React.Fragment key={i}>{children.slice(i, m.index)}</React.Fragment>,
        m[1] === 'doi:' ? <a key={i+m.index} href={`https://doi.org/${m[2]}`}>{m[0]}</a>
        : m[1] === 'http://' ? <a key={i+m.index} href={m[0]}>{m[0]}</a>
        : m[1] === 'https://' ? <a key={i+m.index} href={m[0]}>{m[0]}</a>
        : <React.Fragment key={i+m.index}>{m[0]}</React.Fragment>,
      ]
    }) : { i, el }, { i: 0, el: [] as React.ReactNode[] })
  el.push(<React.Fragment key={i}>{children.slice(i)}</React.Fragment>)
  return el
}
