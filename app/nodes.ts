import * as components from '@/components'
import { MetaNodeAny, MetaNodeData, MetaNodeProcessPrompt, MetaNodeProcessResolve } from '@/spec/metanode'

const nodes = Object.keys(components)
  .map(component_name => components[component_name])
  .filter((component): component is MetaNodeAny => component.spec !== undefined)

export const dataNodes = nodes.reduce<Record<string, MetaNodeData>>((N, component) => component.kind === 'data' ? {...N, [component.spec]: component as any} : N, {})
export const processNodes = nodes.reduce<Record<string, MetaNodeProcessPrompt|MetaNodeProcessResolve>>((N, component) => component.kind === 'process' ? {...N, [component.spec]: component as any} : N, {})
