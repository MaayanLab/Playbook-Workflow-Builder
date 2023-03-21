import React from 'react'
import CGM from 'clustergrammer-gl'
import { v4 as uuidv4 } from 'uuid'

export type Clustergrammer2Data = Omit<Parameters<typeof CGM>[0], 'root'>

export default function ReactCGM(props: Clustergrammer2Data) {
  const ref = React.useRef<HTMLDivElement>(null)
  React.useEffect(() => {
    if (typeof window === 'undefined' || !ref || !ref.current || !props.network) return
    const id = `cgm${uuidv4().replace('-','')}`
    const container = document.createElement('div')
    container.setAttribute('id', id)
    ref.current.appendChild(container)
    const {
      viz_width = 900,
      viz_height = 1035,
      use_hzome = true,
      ...props_
    } = props
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
  }, [ref, props])
  return (
    <div ref={ref} className="flex-grow overflow-scroll relative" />
  )
}
