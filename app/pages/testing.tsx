import React from 'react'
import styles from '@/app/styles/testing.module.css'
import dynamic from 'next/dynamic'
import { MetaNodeDataType } from '@/spec/metanode'
import krg from '@/app/krg'
import * as dict from '@/utils/dict'
import * as array from '@/utils/array'

const JsonEditor = dynamic(() => import('@/app/components/JsonEditor'), { ssr: false })

export default function App() {
  const [data, setData] = React.useState<{
    latest: number,
    current: number,
    selected: Record<number, boolean>,
    nodes: Record<number, { id: number, type: string, data: string, prompt?: string }>,
  }>({
    latest: 0,
    current: 0,
    selected: { [0]: true },
    nodes: { [0]: { id: 0, type: '', data: '' } }
  })
  const [loading, setLoading] = React.useState<boolean>(false)
  const appendData = React.useCallback((add: { type: string, data: string, prompt?: string }) => {
    setData(data => data.nodes[data.current].data ? ({
      ...data,
      latest: data.latest+1,
      current: data.latest+1,
      selected: { [data.latest+1]: true },
      nodes: {
        ...data.nodes,
        [data.latest+1]: {
          id: data.latest+1,
          ...add,
        },
      }
     }) : ({
      ...data,
      selected: { [data.current]: true },
      nodes: {
        ...data.nodes,
        [data.current]: {
          ...data.nodes[data.current],
          ...add,
        },
      },
    }))
  }, [setData])
  const currentNode = data.nodes[data.current]
  const dataNode = krg.getDataNode(currentNode.type)
  let dataNodeView
  if (currentNode.prompt) {
    const promptNode = krg.getPromptNode(currentNode.prompt)
    const inputs: Record<string, unknown> = {}
    if (Object.keys(promptNode.inputs).length > 0) {
      const input0 = Object.keys(promptNode.inputs)[0]
      inputs[input0] = array.ensureOne(promptNode.inputs[input0]).codec.decode(currentNode.data)
    }
    const Prompt = promptNode.prompt
    dataNodeView = <Prompt
      inputs={inputs}
      submit={(output) => {
        setData((data) => ({
          ...data,
          nodes: {
            ...data.nodes,
            [data.current]: {
              ...data.nodes[data.current],
              id: data.current,
              type: promptNode.output.spec,
              data: promptNode.output.codec.encode(output)
            },
          },
        }))
      }}
    />
  } else if (dataNode) {
    try {
      dataNodeView = dataNode.view(dataNode.codec.decode(currentNode.data))
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
        <div className="flex flex-row mt-2 gap-1">
          <div className="flex flex-col-reverse justify-end gap-1">
            {dict.values(data.nodes).filter(item => item.type).map((item) => (
              <button
                key={item.id}
                className={`btn btn-sm ${data.selected[item.id] ? 'btn-primary' : 'btn-secondary'} rounded-md px-2 whitespace-nowrap`}
                onClick={evt => {
                  setData(({ selected: { [item.id]: currentlySelected, ...selected }, ...data }) => ({
                    ...data,
                    current: item.id,
                    selected: {
                      ...(evt.shiftKey ? selected : {}),
                      ...(evt.shiftKey && currentlySelected ? {} : { [item.id]: true }),
                    },
                  }))
                }}
              >{item.type}[{item.id}]</button>
            ))}
          </div>
          <div className="flex flex-col gap-1">
            {krg.getNextProcess(currentNode.type).map(proc =>
              <div key={proc.spec} className="whitespace-nowrap">
                {Object.keys(proc.inputs).length > 0 ? <span> =&gt; </span> : null}
                <button
                  className={[
                    Object.values(proc.inputs)
                      .some((i) => array.ensureOne(i).spec === currentNode.type) ? 'font-bold' : '',
                    'btn btn-sm btn-secondary rounded-sm',
                  ].join(' ')}
                  disabled={!array.all(
                    // all inputs should be satisfiable
                    dict.values(proc.inputs).map((value) => {
                      if (Array.isArray(value)) {
                        return dict.keys(data.selected).filter(id => data.nodes[id].type === value[0].spec).length > 1
                      } else {
                        return dict.keys(data.selected).filter(id => data.nodes[id].type === value.spec).length == 1
                      }
                    })
                  )}
                  onClick={async () => {
                    if ('prompt' in proc) {
                      appendData({
                        type: proc.output.spec,
                        data: '',
                        prompt: proc.spec,
                      })
                    } else {
                      setLoading(() => true)
                      const formData = new FormData()
                      dict.items(proc.inputs).forEach(({ key, value }) => {
                        if (Array.isArray(value)) {
                          dict.keys(data.selected).filter(id => data.nodes[id].type === value[0].spec).forEach(id => {
                            formData.append(key, data.nodes[id].data)
                          })
                        } else {
                          dict.keys(data.selected).filter(id => data.nodes[id].type === value.spec).forEach(id => {
                            formData.append(key, data.nodes[id].data)
                          })
                        }
                      })
                      const req = await fetch(`/api/resolver/${proc.spec}`, {
                        method: 'POST',
                        body: formData,
                      })
                      const res = await req.json()
                      appendData({
                        type: proc.output.spec,
                        data: res,
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
            value={currentNode.data}
            onValueChange={value => {
              setData(data => ({ ...data,
                nodes: {
                  ...data.nodes,
                  [data.current]: {
                    ...data.nodes[data.current],
                    data: value,
                  },
                },
              }))
            }}
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
              setData(_ => ({
                latest: 0,
                current: 0,
                selected: { [0]: true },
                nodes: { [0]: { id: 0, type: '', data: '' } }
              }))
            }}>Reset</button>
          {krg.getDataNodes()
            .filter((node): node is MetaNodeDataType & { meta: { example: unknown } } => 'example' in node.meta && node.meta.example !== undefined)
            .map(node => (
              <button
                key={node.spec}
                className="btn btn-sm btn-secondary rounded-md p-2"
                onClick={() => {
                  appendData({
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
            value={currentNode.type}
            onChange={evt => {
              setData(data => ({ ...data,
                nodes: {
                  ...data.nodes,
                  [data.current]: {
                    ...data.nodes[data.current],
                    type: evt.target.value,
                    prompt: undefined,
                  },
                },
              }))
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
