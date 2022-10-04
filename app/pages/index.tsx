import React from 'react'
import styles from '@/app/styles/index.module.css'
import dynamic from 'next/dynamic'
import { dataNodes, processNodes } from '@/app/nodes'
import { MetaNode, MetaNodeData } from '@/spec/metanode'

const JsonEditor = dynamic(() => import('@/app/components/JsonEditor'))

export default function App() {
  const [data, setData] = React.useState('')
  const [dataType, setDataType] = React.useState('')
  const [prompt, setPrompt] = React.useState(undefined)
  const dataNode = dataNodes[dataType]
  let dataNodeView
  if (prompt) {
    const Prompt = processNodes[prompt].t.prompt
    dataNodeView = <Prompt submit={(output) => {
      setDataType(processNodes[prompt].t.output.t.spec)
      setData(processNodes[prompt].t.output.t.codec.encode(output))
      setPrompt(undefined)
    }}/>
  } else if (dataNode) {
    try {
      dataNodeView = <div>{dataNode.t.view(dataNode.t.codec.decode(data))}</div>
    } catch (e) {
      dataNodeView = <div>Error rendering {dataNode.t.spec}: {e.toString()}</div>
    }
  }
  return (
    <div className={styles.App}>
      <div className={styles.Process}>
        <h2>Apply Process</h2>
        {Object.values(processNodes).map(proc =>
          <div key={proc.t.spec}>
            {Object.keys(proc.t.inputs).length > 0 ?
              <span>{Object.values(proc.t.inputs).map((i: MetaNode<MetaNodeData>) => i.t.spec).join(', ')} =&gt;&nbsp;</span>
              : null}
            <button
              style={{
                fontWeight: Object.values(proc.t.inputs).some((i: MetaNode<MetaNodeData>) => i.t.spec === dataType) ? 'bold' : 'normal',
              }}
              onClick={async () => {
                if ('prompt' in proc.t) {
                  setPrompt(proc.t.spec)
                } else if ('resolve' in proc.t) {
                  const formData = new FormData()
                  for (const i in proc.t.inputs) {
                    formData.append(i, data)
                    // formData[i] = proc.t.inputs[i].t.codec.encode(data)
                  }
                  const req = await fetch(`/api/${proc.t.spec}`, {
                    method: 'POST',
                    body: formData,
                  })
                  const res = await req.json()
                  setPrompt(undefined)
                  setData(res)
                  setDataType(proc.t.output.t.spec)
                }
              }}
            >{proc.t.spec}</button>
            <span>&nbsp; =&gt; {proc.t.output.t.spec}</span>
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
          Load Example:
          {Object.keys(dataNodes)
            .map(spec => dataNodes[spec])
            .filter((node): node is MetaNode<MetaNodeData<unknown, {example: unknown}>> => (node.t.meta as any).example !== undefined)
            .map(node => (
              <button
                key={node.t.spec}
                onClick={() => {
                  setPrompt(undefined)
                  setDataType(node.t.spec)
                  setData(node.t.codec.encode(node.t.meta.example))
                }}
              >{node.t.spec}</button>
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
        >{Object.values(dataNodes).map(dataNode =>
          <option key={dataNode.t.spec} value={dataNode.t.spec}>{dataNode.t.spec}</option>
        )}</select>
        {dataNodeView ? dataNodeView : null}
      </div>
    </div>
  )
}