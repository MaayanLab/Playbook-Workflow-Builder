import React from 'react'
import styles from '@/app/styles/testing.module.css'
import dynamic from 'next/dynamic'
import krg from '@/app/krg'
import * as dict from '@/utils/dict'
import * as array from '@/utils/array'
import classNames from 'classnames'

const Layout = dynamic(() => import('@/app/fragments/playbook/layout'))
const JsonEditor = dynamic(() => import('@/app/components/JsonEditor'), { ssr: false })

export default function App() {
  const [data, setData] = React.useState<{
    // the latest id to be assigned so we can add new ones to the end, always equal to highest node id
    latest: number,
    // the current node in the current data section
    current: number,
    // the nodes we've selected in the apply process section
    selected: Record<number, boolean>,
    // all the nodes, their types, and data
    nodes: Record<number, { id: number, type: string, data: string, prompt?: { type: string, inputs: Record<string, unknown>, data: string } }>,
  }>({
    latest: 0,
    current: 0,
    selected: { [0]: true },
    nodes: { [0]: { id: 0, type: '', data: '' } }
  })
  const [loading, setLoading] = React.useState<boolean>(false)
  /**
   * This is a helper for adding a new node, if the current node doesn't have a type (empty)
   *  we'll just replace that one, otherwise we create a new node, add it, select it, and make it current.
   */
  const appendData = React.useCallback((add: { type: string, data: string, prompt?: { type: string, inputs: Record<string, unknown>, data: string } }) => {
    setData(data => data.nodes[data.current].type ? ({ // our current node has data, add a new one
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
     }) : ({ // our current node has no data, replace it
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
  let promptNodeView
  let dataNodeView
  if (currentNode.prompt) {
    // if our current node has a prompt, show it
    const promptNode = krg.getPromptNode(currentNode.prompt.type)
    const Prompt = promptNode.prompt
    promptNodeView = <Prompt
      data={currentNode.prompt.data ? JSON.parse(currentNode.prompt.data) : undefined}
      output={currentNode.data ? JSON.parse(currentNode.data) : undefined}
      inputs={currentNode.prompt.inputs}
      submit={async (promptData) => {
        if (!currentNode.prompt) return
        setData((data) => {
          return {
            ...data,
            nodes: {
              ...data.nodes,
              [data.current]: {
                ...data.nodes[data.current],
                prompt: {
                  ...data.nodes[data.current].prompt ?? { type: '', data: '', inputs: {} },
                  data: JSON.stringify(promptNode.codec.encode(promptData)),
                },
              },
            },
          }
        })
        setLoading(() => true)
        try {
          const req = await fetch(`/api/resolver`, {
            method: 'POST',
            headers: {
              'Accept': 'application/json',
              'Content-Type': 'application/json',
            },
            body: JSON.stringify({
              spec: promptNode.spec,
              data: promptNode.codec.encode(promptData),
              inputs: currentNode.prompt.inputs,
            }),
          })
          const res = JSON.stringify(await req.json())
          setData(data => ({
            ...data,
            nodes: {
              ...data.nodes,
              [data.current]: {
                ...data.nodes[data.current],
                data: res,
              },
            },
          }))
        } catch (e) {
          appendData({
            type: 'Error',
            data: JSON.stringify((e as Error).toString()),
          })
        } finally {
          setLoading(() => false)
        }
      }}
    />
  }
  if (currentNode.data && dataNode) {
    try {
      dataNodeView = dataNode.view(JSON.parse(currentNode.data))
    } catch (e) {
      dataNodeView = <div>Error rendering {dataNode.meta.label}: {(e as Error).toString()}</div>
    }
  }
  return (
    <Layout>
      <div className={styles.App}>
        <div className={styles.Process}>
          <progress className="progress w-100" value={loading ? undefined : 0}></progress>
          <div className="prose mb-2">
            <h2>Apply Process</h2>
          </div>
          <div className="flex flex-col gap-1">
            {krg.getNextProcess(currentNode.type).map(proc =>
              <div key={proc.spec} className="whitespace-nowrap">
                <button
                  className={[
                    dict.values(proc.inputs)
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
                      const inputs: Record<string, unknown> = {}
                      dict.items(proc.inputs).forEach(({ key, value }) => {
                        if (Array.isArray(value)) {
                          inputs[key] = []
                          dict.keys(data.selected).filter(id => data.nodes[id].type === value[0].spec).forEach(id => {
                            (inputs[key] as string[]).push(JSON.parse(data.nodes[id].data))
                          })
                        } else {
                          dict.keys(data.selected).filter(id => data.nodes[id].type === value.spec).forEach(id => {
                            inputs[key] = JSON.parse(data.nodes[id].data)
                          })
                        }
                      })
                      appendData({
                        type: proc.output.spec,
                        data: '',
                        prompt: {
                          type: proc.spec,
                          inputs,
                          data: '',
                        },
                      })
                    } else {
                      setLoading(() => true)
                      const inputs: Record<string, unknown> = {}
                      dict.items(proc.inputs).forEach(({ key, value }) => {
                        if (Array.isArray(value)) {
                          dict.keys(data.selected).filter(id => data.nodes[id].type === value[0].spec).forEach(id => {
                            inputs[key] = JSON.parse(data.nodes[id].data)
                          })
                        } else {
                          dict.keys(data.selected).filter(id => data.nodes[id].type === value.spec).forEach(id => {
                            inputs[key] = JSON.parse(data.nodes[id].data)
                          })
                        }
                      })
                      try {
                        const req = await fetch(`/api/resolver`, {
                          method: 'POST',
                          headers: {
                            'Accept': 'application/json',
                            'Content-Type': 'application/json',
                          },
                          body: JSON.stringify({
                            spec: proc.spec,
                            inputs,
                          }),
                        })
                        const res = JSON.stringify(await req.json())
                        appendData({
                          type: req.status === 200 ? proc.output.spec : 'Error',
                          data: res,
                        })
                      } catch (e) {
                        appendData({
                          type: 'Error',
                          data: JSON.stringify((e as Error).toString()),
                        })
                      } finally {
                        setLoading(() => false)
                      }
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
          <div style={{
            flex: '1 0 auto',
            height: 0,
            overflow: 'auto',
          }}>
            {currentNode.prompt ? <>
              <span className="prose"><h2>Current Prompt Data</h2></span>
              <JsonEditor
                value={currentNode.prompt.data}
                onValueChange={value => {
                  setData(data => ({ ...data,
                    nodes: {
                      ...data.nodes,
                      [data.current]: {
                        ...data.nodes[data.current],
                        prompt: {
                          ...data.nodes[data.current].prompt ?? { type: '', inputs: {} },
                          data: value,
                        },
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
            </> : null}
            <span className="prose"><h2>Current Output</h2></span>
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
          <div className="flex flex-col gap-2">
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
              <button
                className="btn btn-sm btn-secondary rounded-md p-2"
                onClick={() => {
                  appendData({ type: '', data: '' })
                }}>Start</button>
            </div>
            <div className="flex flex-row flex-wrap gap-1">
              {dict.values(data.nodes).filter(item => item.type).map((item) => (
                <button
                  key={item.id}
                  className={classNames('btn btn-sm rounded-md px-2 whitespace-nowrap', { 'btn-primary': data.selected[item.id], 'btn-secondary': !data.selected[item.id]})}
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
              <option disabled value=""></option>
              {krg.getDataNodes().map(dataNode =>
                <option key={dataNode.spec} value={dataNode.spec}>{dataNode.meta.label}</option>
              )}
            </select>
          </div>
          <div className="m-2 flex-grow flex flex-col gap-2">
            {promptNodeView ? promptNodeView : null}
            {dataNodeView ? dataNodeView : null}
          </div>
        </div>
      </div>
    </Layout>
  )
}
