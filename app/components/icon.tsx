import * as React from 'react'
import type { ReactNode } from 'react'
import { ensureArray } from '@/utils/array'
import type { Icon } from '@/icons'

const base_size = 24
const tooltip_offset_top = 5

export default function Icon({
  icon,
  color: color_,
  size: size_,
  title,
  without_svg,
  container,
  container_color,
}: {
  icon?: Icon,
  color?: string,
  size?: number,
  title?: string | null,
  without_svg?: boolean,
  container?: 'circle' | 'square',
  container_color?: string,
}) {
  const tooltipContainerRef = React.useRef<SVGRectElement>(null)
  const tooltipRef = React.useRef<{ tooltip: HTMLDivElement, listener: (evt: MouseEvent) => void }>(null)
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
    if (tooltipContainerRef.current === null || !title || tooltipRef.current !== null) return
    const svgRect = tooltipContainerRef.current.getBoundingClientRect()
    const left = svgRect.left + window.scrollX
    const right = left + svgRect.width
    const top = svgRect.top + window.scrollY
    const bottom = top + svgRect.height
    const tooltip = document.createElement('div')
    tooltip.setAttribute('data-tip', title as string)
    tooltip.className = 'tooltip tooltip-open'
    tooltip.style.pointerEvents = 'none'
    tooltip.style.position = 'absolute'
    tooltip.style.left = `${left}px`
    tooltip.style.top = `${top - tooltip_offset_top}px`
    tooltip.style.width = `${svgRect.width}px`
    tooltip.style.height = `${svgRect.height + tooltip_offset_top}px`
    tooltip.style.zIndex = '1'
    // since mouseleave/mouseout events are unreliable, we instead watch the mouse
    //  position globally and remove it if it leaves the confines of the rectangle
    // finally -- we write it to a ref so we can remove it if the component is unmounted
    const listener = (evt: MouseEvent) => {
      if ((evt.clientX+window.scrollX < left || evt.clientX+window.scrollX > right)
          || (evt.clientY+window.scrollY < top || evt.clientY+window.scrollY > bottom)) {
        document.removeEventListener('mousemove', listener)
        document.body.removeChild(tooltip)
        Object.assign(tooltipRef, { current: null })
      }
    }
    document.body.append(tooltip)
    Object.assign(tooltipRef, { current: { tooltip, listener } })
    document.addEventListener('mousemove', listener)
  }, [tooltipContainerRef, tooltipRef, title])
  React.useEffect(() => () => {
    if (tooltipRef.current !== null) {
      document.removeEventListener('mousemove', tooltipRef.current.listener)
      document.body.removeChild(tooltipRef.current.tooltip)
    }
  }, [])
  const Svg = without_svg ? (
    ({ children }: { children: ReactNode }) =>
      <g transform={`scale(${size})`}>{children}</g>
  ) : (
    ({ children }: { children: ReactNode }) =>
    <svg
      viewBox={`0 0 ${base_size} ${base_size}`}
      style={{ width: `${size*2}rem`, height: `${size*2}rem`, display: 'inline-block' }}
    >{children}</svg>
  )
  const Container = ({ children }: { children: ReactNode }) =>
    container === 'circle' ? <g transform={`translate(${base_size/2} ${base_size/2})`}>
      <circle r={base_size/2} fill={container_color || '#ddd'} />
      <g transform={`translate(-${base_size/2-0.1*base_size} -${base_size/2-0.1*base_size}) scale(0.8)`}>
        {children}
      </g>
    </g>
    : container === 'square' ? <g>
      <rect width={base_size} height={base_size} fill={container_color || '#ddd'} />
      <g transform={`translate(${0.1*base_size} ${0.1*base_size}) scale(0.8)`}>
        {children}
      </g>
    </g>
    : <>{children}</>
  return (
    <Svg>
      <Container>
        {icons.map((({ path, transform, size }, ind) => {
          return (
            <g key={ind} transform={[
              transformer(ind),
              size && size !== base_size && `scale(${base_size/size})`,
            ].filter((v) => !!v).join(' ')}>
              {transform ? <g transform={transform}><path d={path} style={{ fill: color || 'currentcolor' }} /></g>
                : <path d={path} style={{ fill: color || 'currentcolor' }} />}
            </g>
          )
        }))}
      </Container>
      {title !== null ?
        <rect
          ref={tooltipContainerRef}
          fill="transparent"
          onMouseOver={onMouseOver}
          x={0} width={base_size}
          y={0} height={base_size}
        ></rect>
      : null}
    </Svg>
  )
}