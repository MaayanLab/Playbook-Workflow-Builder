import React from 'react'
import styles from '@/app/styles/testing.module.css'
import dynamic from 'next/dynamic'
import { MetaNodeDataType } from '@/spec/metanode'
import krg from '@/app/krg'
import * as dict from '@/utils/dict'

const JsonEditor = dynamic(() => import('@/app/components/JsonEditor'), { ssr: false })

export default function App() {
  const [prev, setPrev] = React.useState<{ type: string, data: string }[]>([])
  const [loading, setLoading] = React.useState<boolean>(false)
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
  const [prompt, setPrompt] = React.useState<string|undefined>(undefined)
  const dataNode = krg.getDataNode(current.type)
  let dataNodeView
  if (prompt) {
    const promptNode = krg.getPromptNode(prompt)
    const inputs: Record<string, unknown> = {}
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
      dataNodeView = dataNode.view(dataNode.codec.decode(current.data))
    } catch (e) {
      dataNodeView = <div>Error rendering {dataNode.meta.label}: {(e as Error).toString()}</div>
    }
  }
  return (
    <div className={`${styles.App} container mx-auto py-2`}>
      <div className={styles.Process}>
        <progress className="progress w-100" value={loading ? undefined : 0}></progress>
        <div className="prose mb-2">
          <h2>Apply Process</h2>
        </div>
        <div className="flex flex-col gap-1">
          {krg.getNextProcess(current.type).map(proc =>
            <div key={proc.spec} className="whitespace-nowrap">
              {Object.keys(proc.inputs).length > 0 ? (
                <>
                  <span className="btn btn-sm btn-secondary rounded-full">{dict.values(proc.inputs).map((i) => i.meta.label).join(', ')}</span>
                  <span> =&gt; </span>
                </>
              ) : null}
              <button
                className={[
                  Object.values(proc.inputs)
                    .some((i) => i.spec === current.type) ? 'font-bold' : '',
                  'btn btn-sm btn-secondary rounded-sm',
                ].join(' ')}
                onClick={async () => {
                  if ('prompt' in proc) {
                    setPrompt(proc.spec)
                  } else {
                    setLoading(() => true)
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
                    setLoading(() => false)
                  }
                }}
              >{proc.meta.label}</button>
              <span> =&gt; </span>
              <span className="btn btn-sm btn-secondary rounded-full">{proc.output.meta.label}</span>
            </div>
          )}
        </div>
      </div>
      <div className={styles.Data}>
        <div className="prose">
          <h2>Current Data</h2>
        </div>
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
        <div className="flex flex-row flex-wrap mt-1 gap-1">
          <button
            className="btn btn-sm btn-secondary rounded-md p-2"
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
            }}>Previous</button>
          <button
            className="btn btn-sm btn-secondary rounded-md p-2"
            onClick={() => {
              setPrompt(undefined)
              setCurrent({
                type: '',
                data: ''
              })
            }}>Reset</button>
          {krg.getDataNodes()
            .filter((node): node is MetaNodeDataType & { meta: { example: unknown } } => 'example' in node.meta)
            .map(node => (
              <button
                key={node.spec}
                className="btn btn-sm btn-secondary rounded-md p-2"
                onClick={() => {
                  setPrompt(undefined)
                  setCurrent({
                    type: node.spec,
                    data: node.codec.encode(node.meta.example),
                  })
                }}
              >{node.meta.label}</button>
            ))}
        </div>
      </div>
      <div className={styles.View}>
        <div className="prose">
          <h2>Current View</h2>
        </div>
        <div className="form-control w-full max-w-xs">
          <label className="label">
            <span className="label-text">Change the current renderer</span>
          </label>
          <select
            className="select select-bordered"
            value={current.type}
            onChange={evt => {
              setPrompt(undefined)
              setCurrent(({ type: evt.target.value }))
            }}
          >
            <option disabled selected></option>
            {krg.getDataNodes().map(dataNode =>
              <option key={dataNode.spec} value={dataNode.spec}>{dataNode.meta.label}</option>
            )}
          </select>
        </div>
        <div className="m-2 flex-grow flex flex-col">
          {dataNodeView ? dataNodeView : null}
        </div>
      </div>
    </div>
  )
}
