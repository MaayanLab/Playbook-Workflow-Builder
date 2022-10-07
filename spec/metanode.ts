/**
 * Typescript specification for metagraph definitions.
 */

import { Codec } from '@/spec/codec'
import React from 'react'

export type MetaNodeData<C extends {} = unknown, M extends {} = unknown> = {
  spec: string,
  kind: 'data',
  meta: M,
  codec: Codec<C>,
  view(props: C): React.ReactElement,
}

export type MetaNodeProcess<M extends {} = unknown> = {
  spec: string,
  kind: 'process',
  meta: M,
}

export type MetaNodeProcessPrompt<IN = unknown, OUT = unknown, M extends {} = unknown> = MetaNodeProcess<M> & {
  inputs: { [K in keyof IN]: MetaNodeData<IN[K]> },
  output: MetaNodeData<OUT>,
  prompt(props: {
    inputs: {[K in keyof IN]: IN[K]},
    submit: (output: OUT) => void,
  }): React.ReactElement
}

export type MetaNodeProcessResolve<IN = unknown, OUT = unknown, M extends {} = unknown> = MetaNodeProcess<M> & {
  inputs: { [K in keyof IN]: MetaNodeData<IN[K]> },
  output: MetaNodeData<OUT>,
  resolve(props: {
    inputs: {[K in keyof IN]: IN[K]},
  }): OUT
}

export type MetaNodeAny = MetaNodeData | MetaNodeProcess | MetaNodeProcessPrompt | MetaNodeProcessResolve

// Type helpers

/**
 * Extract meta type from a metanode
 */
export type MetaNodeTypeOfMeta<T> = T extends { meta: infer M } ? M : never
/**
 * Extract data type (defined by codec) from a metanode
 */
export type MetaNodeTypeOfCodec<T> = T extends { codec: Codec<infer C> } ? C : never
/**
 * Extract input types
 */
export type MetaNodeTypeOfInputs<T> = T extends { inputs: infer IN } ? { [K in keyof IN]: MetaNodeTypeOfCodec<IN[K]> } : never
/**
 * Extract output type
 */
export type MetaNodeTypeOfOutput<T> = T extends { output: infer OUT } ? MetaNodeTypeOfCodec<OUT> : never

/**
 * This class is used to help build all the attributes of the MetaNodes in a
 *  type-safe manner
 */
export class MetaNode<T = unknown> {
  constructor(public t: T) {}
  /**
   * Begin creating a node
   */
  static createData<S extends string>(spec: S) {
    return new MetaNode({ spec, kind: 'data' as 'data' })
  }
  /**
   * Meta descriptors associated with this node
   */
  meta<M extends {}>(meta: M) {
    return new MetaNode({...this.t, meta})
  }
  /**
   * A codec for the node's underlying data
   */
  codec<C = undefined>(codec: Codec<C> = { encode: JSON.stringify, decode: JSON.parse }) {
    return new MetaNode({ ...this.t, codec })
  }

  /* Data */
  /**
   * `codec` must be defined before view is definable
   */
  view<V extends (props: MetaNodeTypeOfCodec<T>) => React.ReactElement>(view: V) {
    return new MetaNode({ ...this.t, view })
  }

  /* Process */

  /**
   * Begin creating an process
   */
  static createProcess<S extends string>(spec: S) {
    return new MetaNode({ spec, kind: 'process' as 'process' })
  }
  /**
   * The input types of this action, keys mapped to values
   */
  inputs<IN extends { [K in keyof IN]: { codec: Codec<MetaNodeTypeOfCodec<IN[K]>> } }>(inputs: IN = {} as IN) {
    return new MetaNode({ ...this.t, inputs: inputs })
  }
  /**
   * The output type of this action
   */
  output<OUT extends { codec: Codec<MetaNodeTypeOfCodec<OUT>> }>(output: OUT) {
    return new MetaNode({ ...this.t, output })
  }
  /**
   * `input`/`output` must be defined before prompt is definable
   */
  prompt<P extends (props: { inputs: MetaNodeTypeOfInputs<T>, output?: MetaNodeTypeOfOutput<T>, submit: (output: MetaNodeTypeOfOutput<T>) => void }) => React.ReactElement>(prompt: P) {
    return new MetaNode({ ...this.t, prompt })
  }
  /**
   * `input`/`output` must be defined before resolve is definable
   */
  resolve<R extends (props: { inputs: MetaNodeTypeOfInputs<T> }) => Promise<MetaNodeTypeOfOutput<T>>>(resolve: R) {
    return new MetaNode({ ...this.t, resolve })
  }

  /**
   * Finalize the metanode spec
   */
  build() {
    return this.t
  }
}
