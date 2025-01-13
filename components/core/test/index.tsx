import React from 'react'
import { MetaNode } from '@/spec/metanode'
import { z } from 'zod'

export const TestView = MetaNode('TestView')
  .meta({
    label: 'Test',
    description: 'A test view',
    icon: [],
  })
  .codec(z.string())
  .view(data => (
    <div className="alert alert-info">{data}</div>
  ))
  .build()

export const InitialTestProcess = MetaNode('InitialTestProcess')
  .meta({
    label: 'Test',
    description: 'A test process',
    icon: [],
    hidden: true,
  })
  .codec(z.string())
  .inputs({})
  .output(TestView)
  .prompt(props => {
    React.useEffect(() => {
      if (!props.data) props.submit(`${Math.random()}`, true)
    }, [props.data])
    return <></>
  })
  .resolve(async (props) => `${Math.random()}`)
  .story(props => ({}))
  .build()

export const TestProcess = MetaNode('TestProcess')
  .meta({
    label: 'Test',
    description: 'A test process',
    icon: [],
    hidden: true,
  })
  .inputs({ input: TestView })
  .output(TestView)
  .resolve(async (props) => props.inputs.input)
  .story(props => ({}))
  .build()
