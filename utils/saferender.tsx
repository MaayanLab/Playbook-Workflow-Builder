export default function SafeRender<T>({ component: Component, props }: { component: (props: T) => React.ReactNode, props: T }) {
  try {
    return Component(props)
  } catch (e: any) {
    return <div className="alert alert-error">{e.toString()}</div>
  }
}
