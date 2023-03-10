import React from 'react'
import type { FPL } from '@/core/FPPRG'
import type KRG from '@/core/KRG'
import { useSWRImmutableSticky } from '@/utils/use-sticky'
import fetcher from '@/utils/next-rest-fetcher'
import useSWRMap from '@/utils/swr-map'
import * as dict from '@/utils/dict'

export type Metapath = ReturnType<FPL['toJSON']>

/**
 * Retreive output from the API, decode and return
 */
export function useMetapathOutput(krg: KRG, head: Metapath) {
  const { data: rawOutput, isLoading, error } = useSWRImmutableSticky(() => head ? `/api/db/process/${head.process.id}/output` : undefined)
  const processNode = krg.getProcessNode(head.process.type)
  const outputNode = rawOutput ? krg.getDataNode(rawOutput.type) : processNode.output
  const { output, decodeError } = React.useMemo(() => {
      try {
      return { output: rawOutput && outputNode ? outputNode.codec.decode(rawOutput.value) : rawOutput }
    } catch (e: any) {
      console.error(e)
      return { data: { output: undefined, outputNode }, decodeError: process.env.NODE_ENV === 'production' ? 'An unexpected error occurred.' : e.toString() }
    }
  }, [rawOutput, outputNode])
    return { data: { output, outputNode }, error: error || decodeError, isLoading }
}

/**
 * Retreive inputs to this process from outputs of its inputs
 *  We rely on SWR to help us de-duplicate these requests
 */
export function useMetapathInputs(krg: KRG, head: Metapath) {
  const { data: rawInputs, error, isLoading } = useSWRMap<{ type: string, value: any }>(
    dict.values(head.process.inputs).map(({ id }) => `/api/db/process/${id}/output`),
    fetcher
  )
  const processNode = krg.getProcessNode(head.process.type)
  const { inputs, decodeError } = React.useMemo(() => {
    if (!rawInputs) return { inputs: dict.isEmpty(processNode.inputs) ? {} as Record<string, unknown>: undefined }
    try {
      return {
        inputs: dict.init(
          dict.items(processNode.inputs).map(({ key, value }) => {
            if (Array.isArray(value)) {
              return {
                key,
                value: dict.items(head.process.inputs)
                  .filter(({ key: k }) => k.toString().startsWith(`${key}:`))
                  .map(({ value: { id } }) => value[0].codec.decode(rawInputs[`/api/db/process/${id}/output`].value))
              }
            } else {
              return {
                key,
                value: value.codec.decode(rawInputs[`/api/db/process/${head.process.inputs[key].id}/output`].value)
              }
            }
          })
        ),
      }
    } catch (e: any) {
      console.error(e)
      return { decodeError: process.env.NODE_ENV === 'production' ? 'An unexpected error occurred.' : e.toString() }
    }
  }, [rawInputs, processNode])
  return { data: inputs, error: error || decodeError, isLoading }
}
