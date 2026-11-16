import dynamic from 'next/dynamic'
import type { Options as ReactMarkdownOptions } from 'react-markdown'
const Markdown = dynamic<ReactMarkdownOptions>(() => import('react-markdown'), { ssr: false })
export default Markdown