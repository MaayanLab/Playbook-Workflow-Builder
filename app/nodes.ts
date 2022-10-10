import * as components from '@/components'
import { MetaNodeDataType, MetaNodeGenericType, MetaNodePromptType, MetaNodeResolveType } from '@/spec/metanode'

export const dataNodes = Object.keys(components)
  .map(component_name => components[component_name] as MetaNodeGenericType)
  .reduce<Record<string, MetaNodeDataType>>((N, component) => component.kind === 'data' ? { ...N, [component.spec]: component } : N, {})

export const resolveNodes = Object.keys(components)
  .map(component_name => components[component_name] as MetaNodeGenericType)
  .reduce<Record<string, MetaNodeResolveType>>((N, component) => component.kind === 'process' && 'resolve' in component ? { ...N, [component.spec]: component } : N, {})

export const promptNodes = Object.keys(components)
  .map(component_name => components[component_name] as MetaNodeGenericType)
  .reduce<Record<string, MetaNodePromptType>>((N, component) => component.kind === 'process' && 'prompt' in component ? { ...N, [component.spec]: component } : N, {})
