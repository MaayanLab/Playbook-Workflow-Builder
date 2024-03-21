import { ErrorBoundary, FallbackProps } from "react-error-boundary";

function Fallback({ error }: FallbackProps) {
  return <div className="alert alert-error">{error.message}</div>
}

export default function SafeRender<T>({ component: Component, props }: { component: (props: T) => React.ReactNode, props: T & JSX.IntrinsicAttributes }) {
  return <ErrorBoundary FallbackComponent={Fallback}>
    <Component {...props} />
  </ErrorBoundary>
}
