import React from 'react'
import CGM from 'clustergrammer-gl'
import { v4 as uuidv4 } from 'uuid'

export type Clustergrammer2Data = Omit<Parameters<typeof CGM>[0], 'root'>

export default function ReactCGM({ viz_width = 900, viz_height = 1035, ...props }: Clustergrammer2Data) {
  const ref = React.useRef<HTMLDivElement>(null)
  React.useEffect(() => {
    if (typeof window === 'undefined' || !ref || !ref.current || !props.network) return
    const id = `cgm${uuidv4().replace('-','')}`
    const container = document.createElement('div')
    container.setAttribute('id', id)
    ref.current.appendChild(container)
    const { use_hzome = true, ...props_ } = props
    CGM({
      viz_width,
      viz_height,
      use_hzome,
      ...props_,
      container,
    })
    return () => {
      if (ref && ref.current) ref.current.removeChild(container)
    }
  }, [ref, props, viz_width, viz_height])
  return (
    <div ref={ref} style={{ width: `${viz_width}px`, height: `${viz_height + 150}px`, margin: 'auto' }} />
  )
}
