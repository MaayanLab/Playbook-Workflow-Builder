/**
 * Typescript specification for metagraph definitions.
 */

import { Codec } from '@/spec/codec'

export type MetaNodeData<C extends {} = unknown, M extends {} = unknown> = {
  spec: string,
  kind: 'data',
  meta: M,
  codec: Codec<C>,
  view(props: C): any,
}

export type MetaNodeProcess<M extends {} = unknown> = {
  spec: string,
  kind: 'process',
  meta: M,
}

export type MetaNodeProcessPrompt<IN = unknown, OUT = unknown, M extends {} = unknown> = MetaNodeProcess<M> & {
  inputs: { [K in keyof IN]: MetaNode<MetaNodeData<IN[K]>> },
  output: MetaNode<MetaNodeData<OUT>>,
  prompt(props: {
    inputs: {[K in keyof IN]: IN[K]},
    submit: (output: OUT) => void,
  }): any
}

export type MetaNodeProcessResolve<IN = unknown, OUT = unknown, M extends {} = unknown> = MetaNodeProcess<M> & {
  inputs: { [K in keyof IN]: MetaNode<MetaNodeData<IN[K]>> },
  output: MetaNode<MetaNodeData<OUT>>,
  resolve(props: {
    inputs: {[K in keyof IN]: IN[K]},
  }): OUT
}

export type MetaNodeAny = MetaNodeData | MetaNodeProcess | MetaNodeProcessPrompt | MetaNodeProcessResolve

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
  codec<C extends {} = undefined>(codec: Codec<C> = { encode: JSON.stringify, decode: JSON.parse }) {
    return new MetaNode({ ...this.t, codec })
  }

  /* Data */
  /**
   * `codec` must be defined before view is definable
   */
  view<V extends T extends { codec: Codec<infer C> } ? (props: C) => any : never>(view: V) {
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
  inputs<IN extends { [K in keyof IN]: IN[K] }>(inputs: { [K in keyof IN]: MetaNode<IN[K]> } = {} as { [K in keyof IN]: MetaNode<IN[K]> }) {
    return new MetaNode({ ...this.t, inputs: inputs })
  }
  /**
   * The output type of this action
   */
  output<OUT extends {}>(output: MetaNode<OUT>) {
    return new MetaNode({ ...this.t, output })
  }
  /**
   * `input`/`output` must be defined before prompt is definable
   */
  prompt<
    R extends
      T extends {
        inputs: infer IN,
        output: MetaNode<{ codec: Codec<infer OUT> }>
      } ? (props: {
        inputs: {
          [K in keyof IN]: IN[K] extends MetaNode<{ codec: Codec<infer IN_K> }> ? IN_K : never
        },
        output?: OUT,
        submit: (output: OUT) => void,
      }) => React.ReactElement
      : never
  >(prompt: R) {
    return new MetaNode({ ...this.t, prompt })
  }
  /**
   * `input`/`output` must be defined before resolve is definable
   */
  resolve<
    R extends
      T extends {
        inputs: infer IN,
        output: MetaNode<{ codec: Codec<infer OUT> }>
      } ? (props: {
        inputs: {
          [K in keyof IN]: IN[K] extends MetaNode<{ codec: Codec<infer IN_K> }> ? IN_K : never
        }
      }) => Promise<OUT>
      : never
  >(resolve: R) {
    return new MetaNode({ ...this.t, resolve })
  }

  /**
   * Finalize the metanode spec
   */
  build() {
    return this
  }
}
