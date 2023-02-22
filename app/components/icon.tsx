import * as React from 'react'
import type { ReactNode } from 'react'
import { ensureArray } from '@/utils/array'
import type { Icon } from '@/icons'

const base_size = 24
const tooltip_offset_top = 5

export default function Icon({ icon, color: color_, size: size_, title, without_svg }: { icon?: Icon, color?: string, size?: number, title?: string | null, without_svg?: boolean }) {
  const tooltipContainerRef = React.useRef<SVGRectElement>(null)
  const tooltipRef = React.useRef<HTMLDivElement>(null)
  const color = color_ !== undefined ? color_ : 'auto'
  const size = size_ !== undefined ? size_ : 1
  const icons = ensureArray(icon)
  if (title === undefined) title = icons.filter(({ title }) => title).map(({ title }) => title).join(': ')
  const transformer = React.useCallback((ind: number) => {
    if (icons.length === 1) {
      // no change
      return undefined
    } else if (icons.length === 2) {
      // side by side, vertical alignment
      return `translate(${(base_size/2)*ind} ${base_size/4}) scale(0.5, 0.5)`
    } else if (icons.length === 3) {
      // top left, right, & bottom center
      if (ind < 1) {
        return `translate(${base_size/4} 0) scale(0.5, 0.5)`
      } else {
        return `translate(${(base_size/2)*(ind-1)} ${base_size/2}) scale(0.5, 0.5)`
      }
    } else if (icons.length === 4) {
      // 4 quadrants
      if (ind < 2) {
        return `translate(${(base_size/2)*ind} 0) scale(0.5, 0.5)`
      } else {
        return `translate(${(base_size/2)*(ind-2)} ${base_size/2}) scale(0.5, 0.5)`
      }
    } else {
      throw new Error('unhandled')
    }
  }, [icon])
  const onMouseOver = React.useCallback<React.MouseEventHandler<SVGGElement>>((evt) => {
    if (tooltipContainerRef.current === null) return
    let current: HTMLDivElement
    if (tooltipRef.current === null) {
      current = document.createElement('div')
      Object.assign(tooltipRef, { current })
      document.body.append(current)
    } else {
      current = tooltipRef.current
    }
    const svgRect = tooltipContainerRef.current.getBoundingClientRect()
    current.setAttribute('data-tip', title as string)
    current.className = 'tooltip tooltip-open'
    current.style.pointerEvents = 'none'
    current.style.position = 'absolute'
    current.style.left = `${svgRect.left + window.scrollX}px`
    current.style.top = `${svgRect.top + window.scrollY}px`
    current.style.width = `${svgRect.width}px`
    current.style.height = `${svgRect.height}px`
    current.style.zIndex = '1'
  }, [tooltipContainerRef, title])
  const onMouseOut = React.useCallback<React.MouseEventHandler<SVGGElement>>((evt) => {
    if (tooltipRef.current !== null) {
      try {document.body.removeChild(tooltipRef.current)} catch (e) { console.warn(e) }
      Object.assign(tooltipRef, { current: null })
    }
  }, [tooltipContainerRef, title])
  React.useEffect(() => () => {
    if (tooltipRef.current !== null) {
      try {document.body.removeChild(tooltipRef.current)} catch (e) { console.warn(e) }
    }
    }, [])
  const Svg = without_svg ? (
    ({ children }: { children: ReactNode }) =>
      <>{children}</>
  ) : (
    ({ children }: { children: ReactNode }) =>
    <svg
      viewBox={`0 0 ${base_size} ${base_size}`}
      style={{ width: `${size*2}rem`, height: `${size*2}rem`, display: 'inline-block' }}
    >{children}</svg>
  )
  return (
    <Svg>
      {icons.map((({ path, transform, size }, ind) => {
        return (
          <g key={ind} transform={[
            transformer(ind),
            transform,
            size && size !== base_size && `scale(${base_size/size})`,
          ].filter((v) => !!v).join(' ')}>
            <path d={path} style={{ fill: color || 'currentcolor' }} />
          </g>
        )
      }))}
      {title !== null ?
        <rect
          ref={tooltipContainerRef}
          fill="transparent"
          onMouseOver={title !== null ? onMouseOver : undefined}
          onMouseOut={title !== null ? onMouseOut : undefined}
          x={0} width={base_size}
          y={0} height={base_size}
        ></rect>
      : null}
    </Svg>
  )
}