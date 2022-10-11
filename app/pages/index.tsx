import React from 'react'
import styles from '@/app/styles/index.module.css'
import dynamic from 'next/dynamic'
import { MetaNodeDataType } from '@/spec/metanode'
import krg from '@/app/krg'

const JsonEditor = dynamic(() => import('@/app/components/JsonEditor'))

export default function App() {
  const [data, setData] = React.useState('')
  const [dataType, setDataType] = React.useState('')
  const [prompt, setPrompt] = React.useState(undefined)
  const dataNode = krg.getDataNode(dataType)
  let dataNodeView
  if (prompt) {
    const promptNode = krg.getPromptNode(prompt)
    const inputs = {}
    if (Object.keys(promptNode.inputs).length > 0) {
      const input0 = Object.keys(promptNode.inputs)[0]
      inputs[input0] = promptNode.inputs[input0].codec.decode(data)
    }
    const Prompt = promptNode.prompt
    dataNodeView = <Prompt
      inputs={inputs}
      submit={(output) => {
        setDataType(promptNode.output.spec)
        setData(promptNode.output.codec.encode(output))
        setPrompt(undefined)
      }}
    />
  } else if (dataNode) {
    try {
      dataNodeView = <div>{dataNode.view(dataNode.codec.decode(data))}</div>
    } catch (e) {
      dataNodeView = <div>Error rendering {dataNode.spec}: {e.toString()}</div>
    }
  }
  return (
    <div className={styles.App}>
      <div className={styles.Process}>
        <h2>Apply Process</h2>
        {krg.getNextProcess(dataType).map(proc =>
          <div key={proc.spec}>
            {Object.keys(proc.inputs).length > 0 ?
              <span>{Object.values(proc.inputs).map((i) => i.spec).join(', ')} =&gt;&nbsp;</span>
              : null}
            <button
              style={{
                fontWeight: Object.values(proc.inputs).some((i) => i.spec === dataType) ? 'bold' : 'normal',
              }}
              onClick={async () => {
                if ('prompt' in proc) {
                  setPrompt(proc.spec)
                } else {
                  const formData = new FormData()
                  for (const i in proc.inputs) {
                    formData.append(i, data)
                    // formData[i] = proc.inputs[i].codec.encode(data)
                  }
                  const req = await fetch(`/api/resolver/${proc.spec}`, {
                    method: 'POST',
                    body: formData,
                  })
                  const res = await req.json()
                  setPrompt(undefined)
                  setData(res)
                  setDataType(proc.output.spec)
                }
              }}
            >{proc.spec}</button>
            <span>&nbsp; =&gt; {proc.output.spec}</span>
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
            value={data}
            onValueChange={value => setData(value)}
            style={{
              fontFamily: 'monospace',
              fontSize: 12,
              border: '1px solid black',
            }}
          />
        </div>
        <div className={styles.Examples}>
          <button onClick={() => {
            setPrompt(undefined)
            setDataType('')
            setData('')
          }}>Reset</button>
          {krg.getDataNodes()
            .filter((node): node is MetaNodeDataType & { meta: { example: unknown } } => 'example' in node.meta)
            .map(node => (
              <button
                key={node.spec}
                onClick={() => {
                  setPrompt(undefined)
                  setDataType(node.spec)
                  setData(node.codec.encode(node.meta.example))
                }}
              >{node.spec}</button>
            ))}
        </div>
      </div>
      <div className={styles.View}>
        <h2>Current View</h2>
        <select
          value={dataType}
          onChange={evt => {
            setPrompt(undefined)
            setDataType(evt.target.value)
          }}
        >{krg.getDataNodes().map(dataNode =>
          <option key={dataNode.spec} value={dataNode.spec}>{dataNode.spec}</option>
        )}</select>
        {dataNodeView ? dataNodeView : null}
      </div>
    </div>
  )
}