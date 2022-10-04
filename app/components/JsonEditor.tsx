import Editor from 'react-simple-code-editor'
import { highlight, languages } from 'prismjs/components/prism-core'
import 'prismjs/components/prism-json'
import 'prismjs/themes/prism.css'

export default function JsonEditor(props: Omit<React.ComponentProps<typeof Editor>, 'tabSize' | 'insertSpaces' | 'ignoreTabKey' | 'padding' | 'highlight'>) {
  return (
    <Editor
      highlight={value => highlight(value, languages.json)}
      padding={10}
      {...props}
    />
  )
}