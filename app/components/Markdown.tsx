import dynamic from 'next/dynamic'
import type { ReactMarkdownOptions } from 'react-markdown/lib/react-markdown'
const Markdown = dynamic<ReactMarkdownOptions>(() => import('react-markdown') as any, { ssr: false })
export default Markdown