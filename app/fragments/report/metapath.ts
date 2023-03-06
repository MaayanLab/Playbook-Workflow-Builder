import React from 'react'
import type { FPL } from '@/core/FPPRG'
import type KRG from '@/core/KRG'
import type { DataMetaNode } from '@/spec/metanode'
import { useSWRConfig } from 'swr'
import fetcher from '@/utils/next-rest-fetcher'
import useAsyncEffect from 'use-async-effect'
import * as dict from '@/utils/dict'

export type Metapath = ReturnType<FPL['toJSON']>
export type MetapathOutput = {
  outputNode: DataMetaNode,
  output: any,
}

async function errorableFetch<T>(key: string) {
  try {
    return { data: await fetcher<T>(key), error: undefined }
  } catch (e) {
    return { data: undefined, error: e }
  }
}

export type MetapathOutputs = ReturnType<typeof useMetapathOutputs>
export function useMetapathOutputs(krg: KRG, metapath?: Array<Metapath>) {
  const [outputs, setOutputs] = React.useState<Record<string, {
    data?: MetapathOutput,
    isLoading?: boolean,
    error?: unknown
  }>>({})
  const {cache, mutate} = useSWRConfig()
  useAsyncEffect(async (isMounted) => {
    if (!metapath) {
      setOutputs(() => ({}))
      return
    }
    await Promise.all(
      dict.values(dict.init(metapath.map(head => ({ key: head.process.id, value: head.process })))).map(async (head) => {
        const key = `/api/db/process/${head.id}/output`
        let cachedValue = cache.get(key)
        if (cachedValue === undefined || !cachedValue.data) {
          setOutputs(outputs => ({ ...outputs, [head.id]: { isLoading: true } }))
          cachedValue = await errorableFetch(key)
          if (cachedValue.data) {
            await mutate(key, cachedValue.data)
            if (!isMounted()) return
          }
        }
        const { data: rawOutput, error: outputError } = cachedValue
        if (outputError) {
          setOutputs(outputs => ({ ...outputs, [head.id]: { error: outputError } }))
          return
        }
        const processNode = krg.getProcessNode(head.type)
        const outputNode = rawOutput ? krg.getDataNode(rawOutput.type) : processNode.output
        const output = rawOutput && outputNode ? outputNode.codec.decode(rawOutput.value) : rawOutput
        setOutputs(outputs => ({ ...outputs, [head.id]: { data: { outputNode, output } } }))
      })
    )
  }, [krg, metapath])
  return outputs
}

export function useStory(krg: KRG, metapath?: Array<Metapath>, metapathOutputs: MetapathOutputs = {}) {
  return React.useMemo(() => {
    if (!metapath) return
    return metapath
      .map(head => {
        const processNode = krg.getProcessNode(head.process.type)
        if (!processNode.story) return
        const { data: { outputNode = undefined, output = undefined } = {} } = metapathOutputs[head.process.id] || {}
        try {
          return processNode.story({
            inputs: dict.items(head.process.inputs).map(({ key, value }) => ({
              key, value: metapathOutputs[value.id]?.data?.output
            })) as any,
            output,
          })
        } catch (e) {}
      })
      .filter(output => output)
      .join(' ')
  }, [metapath, metapathOutputs])
}
