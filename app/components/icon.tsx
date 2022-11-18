import { ReactNode } from 'react'
import { ensureArray, MaybeArray } from '@/utils/array'

export default function Icon({ icon, color: color_, size: size_, title, without_svg }: { icon?: MaybeArray<{ path: string, transform?: string, title?: string }>, color?: string, size?: number, title?: string, without_svg?: boolean }) {
  const color = color_ !== undefined ? color_ : 'auto'
  const size = size_ !== undefined ? size_ : 1
  const icons = ensureArray(icon)
  const transformer = (ind: number) => {
    if (icons.length === 1) {
      // no change
      return ''
    } else if (icons.length === 2) {
      // side by side, vertical alignment
      return `translate(${12*ind} 6) scale(0.5, 0.5)`
    } else if (icons.length === 3) {
      // top left, right, & bottom center
      if (ind < 1) {
        return `translate(6 0) scale(0.5, 0.5)`
      } else {
        return `translate(${12*(ind-1)} 12) scale(0.5, 0.5)`
      }
    } else if (icons.length === 4) {
      // 4 quadrants
      if (ind < 2) {
        return `translate(${12*ind} 0) scale(0.5, 0.5)`
      } else {
        return `translate(${12*(ind-2)} 12) scale(0.5, 0.5)`
      }
    } else {
      throw new Error('unhandled')
    }
  }
  const Svg = without_svg ? (
    ({ children }: { children: ReactNode }) =>
      <g>{children}</g>
  ) : (
    ({ children }: { children: ReactNode }) =>
    <svg
      viewBox="0 0 24 24"
      style={{ width: `${size*2}rem`, height: `${size*2}rem` }}
    >{children}</svg>
  )
  if (title === undefined) title = icons.filter(({ title }) => title).map(({ title }) => title).join(': ')
  return (
    <Svg>
      {icons.map((({ path, transform }, ind) => {
        return (
          <g key={ind} transform={transformer(ind)}>
            <g transform={transform}>
              <path d={path} style={{ fill: color || 'currentcolor' }} />
            </g>
          </g>
        )
      }))}
      {title ? <title>{title}</title> : null}
    </Svg>
  )
}