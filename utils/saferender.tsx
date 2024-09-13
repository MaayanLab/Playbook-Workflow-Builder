import dynamic from 'next/dynamic';

const ErrorBoundary = dynamic(() => import('react-error-boundary').then(({ ErrorBoundary }) => ErrorBoundary), { ssr: false })

function fallbackRender({ error }: { error: { message: string } }) {
  return (
    <div className="alert alert-error">Client side error: {error.message}</div>
  )
}

export default function SafeRender<T>({ component: Component, props }: { component: (props: T) => React.ReactNode, props: T }) {
  try {
    return <ErrorBoundary fallbackRender={fallbackRender}>{Component(props)}</ErrorBoundary>
  } catch (e: any) {
    return <div className="alert alert-error">Server side error: {e.toString()}</div>
  }
}
