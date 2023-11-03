import React from 'react'
import { Map } from 'immutable'

const WaypointsContext = React.createContext({
  waypoints: Map<string, { ref: HTMLDivElement, active?: boolean }>(),
  set: (id: string, ref: HTMLDivElement) => {},
  del: (id: string) => {},
})

export function useWaypoints() {
  const { waypoints } = React.useContext(WaypointsContext)
  return waypoints
}

export function Waypoints(props: React.PropsWithChildren<{ className?: string }>) {
  const ref = React.useRef<HTMLDivElement>(null)
  const [waypoints, setWaypoints] = React.useState(Map<string, { ref: HTMLDivElement, active?: boolean }>())
  React.useEffect(() => {
    if (!ref.current) return
    const listener = () => {
      if (!ref.current) return
      const scrollTop = ref.current.scrollTop + window.innerHeight/4
      const scrollBottom = scrollTop + window.innerHeight/2
      setWaypoints(waypoints => waypoints.mapEntries(([key, value]) => {
        const elTop = value.ref.offsetTop
        const elBottom = value.ref.offsetTop + value.ref.clientHeight
        const active = (scrollTop < elBottom && scrollBottom > elTop)
        return [key, {
          ref: value.ref,
          active,
        }]
      }))
    }
    ref.current.addEventListener('scroll', listener)
    return () => {
      ref.current?.removeEventListener('scroll', listener)
    }
  }, [ref])
  return (
    <div ref={ref} className={props.className}>
      <WaypointsContext.Provider value={{
        waypoints,
        set: (id: string, ref: HTMLDivElement) => {
          setWaypoints(waypoints => waypoints.set(id, { ref }))
        },
        del: (id: string) => {
          setWaypoints(waypoints => waypoints.delete(id))
        },
      }}>
        {props.children}
      </WaypointsContext.Provider>
    </div>
  )
}

export function Waypoint({ id, children, ...props }: React.PropsWithChildren<{ id: string } & React.HTMLAttributes<HTMLDivElement>>) {
  const ref = React.useRef<HTMLDivElement>(null)
  const waypoints = React.useContext(WaypointsContext)
  React.useEffect(() => {
    if (ref.current === null) return
    waypoints.set(id, ref.current)
    return () => {waypoints.del(id)}
  }, [id, ref])
  return <div ref={ref} {...props}>{children}</div>
}
