import React from 'react'
import styles from '@/app/styles/testing.module.css'
import dynamic from 'next/dynamic'
import { MetaNodeDataType } from '@/spec/metanode'
import krg from '@/app/krg'
import * as dict from '@/utils/dict'

const JsonEditor = dynamic(() => import('@/app/components/JsonEditor'), { ssr: false })
const Button = dynamic(() => import('@blueprintjs/core').then(({ Button }) => Button))

export default function App() {
  const [prev, setPrev] = React.useState([])
  const [current, setCurrent_] = React.useState({ type: '', data: '' })
  const setCurrent = React.useCallback((current: { type?: string, data?: string }) => {
    setCurrent_(current_ => {
      setPrev(prev => [...prev, current_])
      return {
        type: current.type === undefined ? current_.type : current.type,
        data: current.data === undefined? current_.data : current.data,
      }
    })
  }, [setPrev, setCurrent_])
  const [prompt, setPrompt] = React.useState(undefined)
  const dataNode = krg.getDataNode(current.type)
  let dataNodeView
  if (prompt) {
    const promptNode = krg.getPromptNode(prompt)
    const inputs = {}
    if (Object.keys(promptNode.inputs).length > 0) {
      const input0 = Object.keys(promptNode.inputs)[0]
      inputs[input0] = promptNode.inputs[input0].codec.decode(current.data)
    }
    const Prompt = promptNode.prompt
    dataNodeView = <Prompt
      inputs={inputs}
      submit={(output) => {
        setCurrent({
          type: promptNode.output.spec,
          data: promptNode.output.codec.encode(output)
        })
        setPrompt(undefined)
      }}
    />
  } else if (dataNode) {
    try {
      dataNodeView = <div>{dataNode.view(dataNode.codec.decode(current.data))}</div>
    } catch (e) {
      dataNodeView = <div>Error rendering {dataNode.meta.label}: {e.toString()}</div>
    }
  }
  return (
    <div className={styles.App}>
      <div className={styles.Process}>
        <h2>Apply Process</h2>
        {krg.getNextProcess(current.type).map(proc =>
          <div key={proc.spec}>
            {Object.keys(proc.inputs).length > 0 ? (
              <>
                <span className="bg-secondary rounded-full p-3">{dict.values(proc.inputs).map((i) => i.meta.label).join(', ')}</span>
                <span> =&gt; </span>
              </>
            ) : null}
            <Button
              className={[
                Object.values(proc.inputs)
                  .some((i) => i.spec === current.type) ? 'font-bold' : '',
                'bg-primary rounded-sm p-2',
              ].join(' ')}
              onClick={async () => {
                if ('prompt' in proc) {
                  setPrompt(proc.spec)
                } else {
                  const formData = new FormData()
                  for (const i in proc.inputs) {
                    formData.append(i, current.data)
                    // formData[i] = proc.inputs[i].codec.encode(data)
                  }
                  const req = await fetch(`/api/resolver/${proc.spec}`, {
                    method: 'POST',
                    body: formData,
                  })
                  const res = await req.json()
                  setPrompt(undefined)
                  setCurrent({
                    type: proc.output.spec,
                    data: res
                  })
                }
              }}
            >{proc.meta.label}</Button>
            <span> =&gt; </span>
            <span className="bg-secondary rounded-full p-3">{proc.output.meta.label}</span>
          </div>
        )}
      </div>
      <div className={styles.Data}>
        <h2>Current Data</h2>
        <div style={{
          flex: '1 0 auto',
          height: 0,
          overflow: 'auto',
        }}>
          <JsonEditor
            value={current.data}
            onValueChange={value => setCurrent({ data: value })}
            style={{
              fontFamily: 'monospace',
              fontSize: 12,
              border: '1px solid black',
            }}
          />
        </div>
        <div className={styles.Examples}>
          <Button
            className="bg-primary rounded-md p-2"
            onClick={() => {
              setPrompt(undefined)
              setPrev(prev => {
                const _prev = [...prev]
                setCurrent_(_prev.pop() || {
                  type: '',
                  data: ''
                })
                return _prev
            })
            }}>Previous</Button>
          <Button
            className="bg-primary rounded-md p-2"
            onClick={() => {
              setPrompt(undefined)
              setCurrent({
                type: '',
                data: ''
              })
            }}>Reset</Button>
          {krg.getDataNodes()
            .filter((node): node is MetaNodeDataType & { meta: { example: unknown } } => 'example' in node.meta)
            .map(node => (
              <Button
                key={node.spec}
                className="bg-primary rounded-md p-2"
                onClick={() => {
                  setPrompt(undefined)
                  setCurrent({
                    type: node.spec,
                    data: node.codec.encode(node.meta.example),
                  })
                }}
              >{node.meta.label}</Button>
            ))}
        </div>
      </div>
      <div className={styles.View}>
        <h2>Current View</h2>
        <select
          value={current.type}
          onChange={evt => {
            setPrompt(undefined)
            setCurrent(({ type: evt.target.value }))
          }}
        >{krg.getDataNodes().map(dataNode =>
          <option key={dataNode.spec} value={dataNode.spec}>{dataNode.meta.label}</option>
        )}</select>
        {dataNodeView ? dataNodeView : null}
      </div>
    </div>
  )
}