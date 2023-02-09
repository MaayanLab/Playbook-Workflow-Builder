import * as React from 'react'
import type { ReactNode } from 'react'
import { ensureArray } from '@/utils/array'
import type { Icon } from '@/icons'

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
      document.body.removeChild(tooltipRef.current)
      Object.assign(tooltipRef, { current: null })
    }
  }, [tooltipContainerRef, title])
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
      {title !== null ?
        <rect
          ref={tooltipContainerRef}
          fill="transparent"
          onMouseOver={title !== null ? onMouseOver : undefined}
          onMouseOut={title !== null ? onMouseOut : undefined}
          x={-4}
          y={-4}
          style={{ width: `${size*2}rem`, height: `${size*2}rem` }}
        ></rect>
      : null}
    </Svg>
  )
}