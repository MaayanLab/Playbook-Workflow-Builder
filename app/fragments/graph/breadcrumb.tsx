import React from 'react'
import { BreadcrumbNode, useBreadcrumbPosition } from "@/app/fragments/breadcrumbs"
import Icon from '@/app/components/icon'
import { type Icon as IconT } from '@/icons'
import * as array from '@/utils/array'
import { Metapath, useResolved } from "@/app/fragments/metapath"

type GraphPosition = {
  head?: Metapath,
  active: boolean,
}
type BreadcrumbStyle = {
  label: string,
  icon: IconT,
}

export function DataBreadcrumb(props: GraphPosition & BreadcrumbStyle & BreadcrumbNode) {
  const position = useBreadcrumbPosition(props)
  const { data } = useResolved(props.head?.process.id)
  const title = React.useMemo(() =>
    props.label || array.ensureArray(props.icon).map(({ title }) => title).join(': ') || props.id,
    [props.label, props.icon, props.id])
  if (!position) return null
  return (
    <g
      className="cursor-pointer"
      transform={`translate(${position.x*1.25} ${position.y*1.25})`}
      onClick={props.onClick}
    >
      <circle
        fill={props.active ? '#B3CFFF'
        : data?.type === 'error' ? '#ff0000'
        : 'lightgrey'}
        r={0.5}
        cx={0}
        cy={0}
      />
      <g transform={`scale(0.035 0.035) translate(-12 -12)`}>
        <Icon icon={props.icon} title={title} without_svg />
      </g>
    </g>
  )
}

export function ProcessBreadcrumb(props: GraphPosition & BreadcrumbStyle & BreadcrumbNode) {
  const position = useBreadcrumbPosition(props)
  const { data } = useResolved(props.head?.process.id)
  const title = React.useMemo(() =>
    props.label || array.ensureArray(props.icon).map(({ title }) => title).join(': ') || props.id,
    [props.label, props.icon, props.id])
  if (!position) return null
  return (
    <g
      className="cursor-pointer"
      transform={`translate(${position.x*1.25} ${position.y*1.25})`}
      onClick={props.onClick}
    >
      <rect
        fill={props.active ? '#B3CFFF'
        : data?.type === 'resolving' ? '#88ff88'
        : data?.type === 'resolved' && data.data.data === null ? '#ffffbb'
        : data?.type === 'resolved' ? '#88ff88'
        : data?.type === 'error' ? '#ffbbbb'
        : 'lightgrey'}
        x={-0.5}
        y={-0.5}
        width={1}
        height={1}
      />
      {data?.type === 'resolving' ?
        <rect
          fill={'#ddffcc'}
          x={-0.5 + data.percent}
          y={-0.5}
          width={1.0 - data.percent}
          height={1}
        /> : null}
      <g transform={`scale(0.035 0.035) translate(-12 -12)`}>
        <Icon icon={props.icon} title={title} without_svg />
      </g>
    </g>
  )
}
