/**
 * Typescript specification for metagraph definitions.
 */

import { Codec } from '@/spec/codec'
import React from 'react'
import { z } from 'zod'
import type { Icon } from '@/icons'
import codecFrom from '@/utils/zod-codec'
import type { MaybeArray, ValuesOfMaybeArray } from '@/utils/types'

/**
 * These type helpers beak up the meta node definitions in such a way that they can be
 *  built up in a type safe fashion.
 */

export type MetaNodeMetadata<T = unknown> = {
  label: string,
  description: string,
  icon?: Icon,
  color?: string,
  default?: T,
  example?: T,
  pagerank?: number
  tags?: Record<string, Record<string, number>>,
}
export type MetaNodeGeneric<T = unknown> = { kind: 'data', spec: string, meta: MetaNodeMetadata<T>, codec: Codec<T> }
export type MetaNodeExtractKind<T = MetaNodeGeneric> = T extends MaybeArray<{ kind: infer Kind }> ? Kind : never
export type MetaNodeExtractSpec<T = MetaNodeGeneric> = T extends MaybeArray<{ spec: infer Spec }> ? Spec : never
export type MetaNodeExtractMeta<T = MetaNodeGeneric> = T extends MaybeArray<{ meta: infer Meta }> ? Meta : never
export type MetaNodeExtractData<T = MetaNodeGeneric> =
  T extends Array<{ codec: Codec<infer Data> }> ? Data[]
  : T extends { codec: Codec<infer Data> } ? Data
  : never
export type MetaNodeView<T = MetaNodeGenericData> = (props: MetaNodeExtractData<T>) => React.ReactElement
export type MetaNodeWithKind<T = MetaNodeGenericData> = { kind: MetaNodeExtractKind<T> }
export type MetaNodeWithSpec<T = MetaNodeGenericData> = { spec: MetaNodeExtractSpec<T> }
export type MetaNodeWithMeta<T = MetaNodeGenericData> = { meta: MetaNodeExtractMeta<T> }
export type MetaNodeWithCodec<T = MetaNodeGenericData> = { codec: Codec<MetaNodeExtractData<T>> }
export type MetaNodeWithView<T = MetaNodeGenericData> = T extends { view: MetaNodeView<T> } ? T : never
export type MetaNodeGenericData = MetaNodeGeneric & { view: MetaNodeView<MetaNodeGeneric> }
export type MetaNodeDataType<T = MetaNodeGenericData> = MetaNodeWithKind<T> & MetaNodeWithSpec<T> & MetaNodeWithMeta<T> & MetaNodeWithCodec<T> & MetaNodeWithView<T>
export type MetaNodeGenericProcess = { kind: 'process', spec: string, meta: MetaNodeMetadata, inputs: Record<string, MaybeArray<MetaNodeGenericData>>, output: MetaNodeGenericData }
export type MetaNodeInputs<T = Record<string, MaybeArray<MetaNodeGenericData>>> = T extends { [K in keyof T]: MaybeArray<MetaNodeDataType<ValuesOfMaybeArray<T[K]>>> } ? T : never
export type MetaNodeExtractInputs<T = MetaNodeGenericProcess> = T extends { inputs: MetaNodeInputs<infer Inputs> } ? Inputs : never
export type MetaNodeOutput<T = MetaNodeGenericData> = T extends MetaNodeDataType<T> ? T : never
export type MetaNodeExtractOutput<T = MetaNodeGenericProcess> = T extends { output: MetaNodeOutput<infer Output> } ? Output : never
export type MetaNodeWithInputs<T = MetaNodeGenericProcess> = T extends { inputs: MetaNodeExtractInputs<T> } ? T : never
export type MetaNodeWithOutput<T = MetaNodeGenericProcess> = T extends { output: MetaNodeExtractOutput<T> } ? T : never
export type MetaNodeExtractInputsData<T = MetaNodeGenericProcess> = { [K in keyof MetaNodeExtractInputs<T>]: MetaNodeExtractData<MetaNodeExtractInputs<T>[K]> }
export type MetaNodeExtractOutputData<T = MetaNodeGenericProcess> = MetaNodeExtractData<MetaNodeExtractOutput<T>>
export type MetaNodePrompt<T = MetaNodeGenericProcess> = (props: { inputs: MetaNodeExtractInputsData<T>, output?: MetaNodeExtractOutputData<T>, submit: (output: MetaNodeExtractOutputData<T>) => void }) => React.ReactElement
export type MetaNodeResolve<T = MetaNodeGenericProcess> = (props: { inputs: MetaNodeExtractInputsData<T> }) => Promise<MetaNodeExtractOutputData<T>>
export type MetaNodeWithPrompt<T = MetaNodeGenericProcess> = T extends { prompt: MetaNodePrompt<T> } ? T : never
export type MetaNodeWithResolve<T = MetaNodeGenericProcess> = T extends { resolve: MetaNodeResolve<T> } ? T : never
export type MetaNodeGenericPrompt = MetaNodeGenericProcess & { prompt: MetaNodePrompt<MetaNodeGenericProcess> }
export type MetaNodePromptType<T = MetaNodeGenericPrompt> = MetaNodeWithKind<T> & MetaNodeWithSpec<T> & MetaNodeWithMeta<T> & MetaNodeWithInputs<T> & MetaNodeWithOutput<T> & MetaNodeWithPrompt<T>
export type MetaNodeGenericResolve = MetaNodeGenericProcess & { resolve: MetaNodeResolve<MetaNodeGenericProcess> }
export type MetaNodeResolveType<T = MetaNodeGenericResolve> = MetaNodeWithKind<T> & MetaNodeWithSpec<T> & MetaNodeWithMeta<T> & MetaNodeWithInputs<T> & MetaNodeWithOutput<T> & MetaNodeWithResolve<T>
export type MetaNodeGenericType = MetaNodeDataType | MetaNodePromptType | MetaNodeResolveType
export type MetaNodeType<T> = MetaNodeDataType<T> | MetaNodePromptType<T> | MetaNodeResolveType<T>

/**
 * This class is used to help build all the attributes of the MetaNodes in a
 *  type-safe manner
 */
export class MetaNode<T = unknown> {
  constructor(public t: T) { }
  /**
   * Begin creating a node
   */
  static createData(spec: string) {
    return new MetaNode({ spec, kind: 'data' as 'data' })
  }
  /**
   * Meta descriptors associated with this node
   */
  meta<M extends MetaNodeMetadata>(meta: M & MetaNodeMetadata) {
    return new MetaNode({ ...this.t, meta })
  }
  /**
   * A codec for the node's underlying data
   */
  codec<C = undefined>(codec: z.ZodType<C> | Codec<C> = { encode: JSON.stringify, decode: JSON.parse }) {
    const codec_: Codec<C> = ('parse' in codec) ? codecFrom(codec) : codec
    return new MetaNode({ ...this.t, codec: codec_ })
  }

  /* Data */
  /**
   * `codec` must be defined before view is definable
   */
  view(view: MetaNodeView<T>) {
    return new MetaNode({ ...this.t, view })
  }

  /* Process */

  /**
   * Begin creating an process
   */
  static createProcess(spec: string) {
    return new MetaNode({ spec, kind: 'process' as 'process' })
  }
  /**
   * The input types of this action, keys mapped to values
   */
  inputs<IN = {}>(inputs: MetaNodeInputs<IN> = {} as MetaNodeInputs<IN>) {
    return new MetaNode({ ...this.t, inputs: inputs })
  }
  /**
   * The output type of this action
   */
  output<OUT>(output: MetaNodeOutput<OUT>) {
    return new MetaNode({ ...this.t, output })
  }
  /**
   * `input`/`output` must be defined before prompt is definable
   */
  prompt(prompt: MetaNodePrompt<T>) {
    return new MetaNode({ ...this.t, prompt })
  }
  /**
   * `input`/`output` must be defined before resolve is definable
   */
  resolve(resolve: MetaNodeResolve<T>) {
    return new MetaNode({ ...this.t, resolve })
  }

  /**
   * Finalize the metanode spec
   */
  build(this: { t: MetaNodeType<T> }) {
    return this.t
  }
}
