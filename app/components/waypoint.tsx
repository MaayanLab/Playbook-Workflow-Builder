/**
 * Reference containers throughout the page which we can:
 *  - persistently scroll to with an anchor
 *  - detect when scrolled into view
 * 
 * Some special waypoints are also assumed:
 * - top: the top of the page
 * - head: the bottom of any sticky header (we calculate & inject the necessary scrollMarginTop)
 * - bottom: bottom of the page
 */
import React from 'react'
import { Map } from 'immutable'

const WaypointsContext = React.createContext({
  refs: Map<string, HTMLDivElement>(),
  waypoints: Map<string, { scrollMarginTop: string, active: boolean }>(),
  set: (id: string, ref: HTMLDivElement) => {},
  del: (id: string) => {},
})

export function useWaypoints() {
  const { refs, waypoints } = React.useContext(WaypointsContext)
  const scrollTo = React.useCallback((id: string) => {
    history.pushState(null, document.title, id === 'top' ? ' ' : `#${id}`)
    window.dispatchEvent(new Event('hashchange'))
  }, [])
  return {
    waypoints: waypoints.mapEntries(([key, value]) => [key, { ...value, ref: refs.get(key) }]),
    scrollTo,
  }
}

export function Waypoints(props: React.PropsWithChildren<{ className?: string }>) {
  const ref = React.useRef<HTMLDivElement>(null)
  const first = React.useRef(true)
  const [waypointRefs, setWaypointRefs] = React.useState(Map<string, HTMLDivElement>())
  const [waypoints, setWaypoints] = React.useState(Map<string, { scrollMarginTop: string, active: boolean }>())
  const set = React.useCallback((id: string, setRef: HTMLDivElement) => {
    setWaypointRefs(waypoints => waypoints.set(id, setRef))
  }, [])
  const del = React.useCallback((id: string) => {
    setWaypointRefs(waypoints => waypoints.delete(id))
  }, [])
  React.useEffect(() => {
    if (!ref.current || waypointRefs.size === 0) return
    const listener = () => {
      if (!ref.current) return
      const scrollTop = ref.current.scrollTop + window.innerHeight/4
      const scrollBottom = scrollTop + window.innerHeight/2
      const head = waypointRefs.get('head')
      const scrollMarginTop = head ? `${head.clientHeight + 2}px` : `0px`
      setWaypoints(() => waypointRefs.mapEntries(([key, value]) => {
        const elTop = value.offsetTop
        const elBottom = value.offsetTop + value.clientHeight
        const active = (scrollTop < elBottom && scrollBottom > elTop)
        return [key, { scrollMarginTop, active }]
      }))
    }
    ref.current.addEventListener('scroll', listener)
    listener()
    return () => {
      ref.current?.removeEventListener('scroll', listener)
    }
  }, [ref, waypointRefs])
  React.useEffect(() => {
    if (!ref.current || waypointRefs.size === 0) return
    const listener = () => {
      const hash = (window.location.hash || '#').slice(1) || 'top'
      const waypointRef = waypointRefs.get(hash)
      if (waypointRef) {
        setTimeout(() => {
          first.current = false
          waypointRef.scrollIntoView({ behavior: 'smooth', block: 'start', inline: 'nearest' })
        }, 200)
      }
    }
    window.addEventListener('hashchange', listener)
    if (first.current) listener()
    return () => {
      window.removeEventListener('hashchange', listener)
    }
  }, [ref, first, waypointRefs])
  return (
    <div ref={ref} className={props.className}>
      <WaypointsContext.Provider value={{ refs: waypointRefs, waypoints, set, del }}>
        {props.children}
      </WaypointsContext.Provider>
    </div>
  )
}

export function Waypoint({ id, children, ...props }: React.PropsWithChildren<{ id: string } & React.HTMLAttributes<HTMLDivElement>>) {
  const ref = React.useRef<HTMLDivElement>(null)
  const { waypoints, set, del } = React.useContext(WaypointsContext)
  const waypoint = React.useMemo(() => waypoints.get(id), [id, waypoints])
  React.useEffect(() => {
    if (ref.current === null) return
    set(id, ref.current)
    return () => { del(id)}
  }, [id, set, del, ref])
  return <div
    ref={ref}
    {...{
      ...props,
      style: {
        ...props.style,
        scrollMarginTop: props.style?.scrollMarginTop ?? waypoint?.scrollMarginTop
      }
    }}
  >{children}</div>
}
