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
        const cachedValue = cache.get(key)
        if (cachedValue !== undefined) {
          setOutputs(outputs => ({ ...outputs, [head.id]: { data: cachedValue as any } }))
          return
        }
        setOutputs(outputs => ({ ...outputs, [head.id]: { isLoading: true } }))
        const { data: rawOutput, error: outputError } = await errorableFetch<any>(key)
        if (!isMounted()) return
        if (outputError) {
          setOutputs(outputs => ({ ...outputs, [head.id]: { error: outputError } }))
          return
        }
        const processNode = krg.getProcessNode(head.type)
        const outputNode = rawOutput ? krg.getDataNode(rawOutput.type) : processNode.output
        const output = rawOutput && outputNode ? outputNode.codec.decode(rawOutput.value) : rawOutput
        await mutate(key, { outputNode, output })
        if (!isMounted()) return
        setOutputs(outputs => ({ ...outputs, [head.id]: { data: { outputNode, output } } }))
      })
    )
  }, [krg, metapath])
  return outputs
}
