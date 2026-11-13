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
import { mapDelete, mapSet, mapValues } from '@/utils/map'

type WaypointState = { scrollMarginTop: string, active: boolean }

const WaypointsContext = React.createContext({
  refs: new Map() as ReadonlyMap<string, HTMLDivElement>,
  waypoints: new Map() as ReadonlyMap<string, WaypointState>,
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
    waypoints: React.useMemo(() => mapValues(waypoints, (value, key) => ({ ...value, ref: refs.get(key) })), [refs, waypoints]),
    scrollTo,
  }
}

export function Waypoints(props: React.PropsWithChildren<{ className?: string }>) {
  const ref = React.useRef<HTMLDivElement>(null)
  const first = React.useRef(true)
  const [waypointRefs, setWaypointRefs] = React.useState<ReadonlyMap<string, HTMLDivElement>>(() => new Map())
  const [waypoints, setWaypoints] = React.useState<ReadonlyMap<string, WaypointState>>(() => new Map())
  const set = React.useCallback((id: string, setRef: HTMLDivElement) => {
    setWaypointRefs(waypoints => mapSet(waypoints, id, setRef))
  }, [])
  const del = React.useCallback((id: string) => {
    setWaypointRefs(waypoints => mapDelete(waypoints, id))
  }, [])
  React.useEffect(() => {
    if (!ref.current || waypointRefs.size === 0) return
    // the page scrolls within a descendant (scroll events don't bubble, but they can be captured),
    //  waypoints are active when they're in the middle half of the window
    const listener = () => {
      const scrollTop = window.innerHeight/4
      const scrollBottom = scrollTop + window.innerHeight/2
      const head = waypointRefs.get('head')
      const scrollMarginTop = head ? `${head.clientHeight + 2}px` : `0px`
      const next = mapValues(waypointRefs, (value) => {
        const { top: elTop, bottom: elBottom } = value.getBoundingClientRect()
        const active = (scrollTop < elBottom && scrollBottom > elTop)
        return { scrollMarginTop, active }
      })
      // keep the state when nothing changed (this runs on every scroll event)
      setWaypoints(prev => prev.size === next.size && Array.from(next).every(([key, value]) => {
        const prevValue = prev.get(key)
        return prevValue?.active === value.active && prevValue.scrollMarginTop === value.scrollMarginTop
      }) ? prev : next)
    }
    const el = ref.current
    el.addEventListener('scroll', listener, { capture: true, passive: true })
    listener()
    return () => {
      el.removeEventListener('scroll', listener, { capture: true })
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
