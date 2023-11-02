import { BreadcrumbNode, useBreadcrumbPosition } from "@/app/fragments/breadcrumbs"
import Icon from '@/app/components/icon'
import { type Icon as IconT } from '@/icons'
import * as array from '@/utils/array'

export type BreadcrumbStyle = {
  kind: string,
  label: string,
  color: string,
  icon: IconT,
}

export function Breadcrumb(props: BreadcrumbStyle & BreadcrumbNode) {
  const position = useBreadcrumbPosition(props)
  if (!position) return null
  const title = props.label || array.ensureArray(props.icon).map(({ title }) => title).join(': ') || props.id
  return (
    <g
      className="cursor-pointer"
      transform={`translate(${position.x*1.25} ${position.y*1.25})`}
      onClick={props.onClick}
    >
      {props.kind === 'data' ? (
        <circle
          fill={props.color}
          r={0.5}
          cx={0}
          cy={0}
        />
      ) : (
        <rect
          fill={props.color}
          x={-0.5}
          y={-0.5}
          width={1}
          height={1}
        />
      )}
      <g transform={`scale(0.035 0.035) translate(-12 -12)`}>
        <Icon icon={props.icon} title={title} without_svg />
      </g>
    </g>
  )
}
