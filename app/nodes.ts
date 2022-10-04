import * as components from '@/components'
import { MetaNode, MetaNodeAny, MetaNodeData, MetaNodeProcessPrompt, MetaNodeProcessRun } from '@/spec/metanode'

const nodes = Object.keys(components)
  .map(component_name => components[component_name])
  .filter((component): component is MetaNode<MetaNodeAny> => component instanceof MetaNode)

export const dataNodes = nodes.reduce<Record<string, MetaNode<MetaNodeData>>>((N, component) => component.t.kind === 'data' ? {...N, [component.t.spec]: component as any} : N, {})
export const processNodes = nodes.reduce<Record<string, MetaNode<MetaNodeProcessPrompt | MetaNodeProcessRun>>>((N, component) => component.t.kind === 'process' ? {...N, [component.t.spec]: component as any} : N, {})
