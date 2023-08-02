import React from 'react'
import type { FPL } from '@/core/FPPRG'
import type KRG from '@/core/KRG'
import useSWRImmutable from 'swr/immutable'
import fetcher from '@/utils/next-rest-fetcher'
import useSWRMap from '@/utils/swr-map'
import * as dict from '@/utils/dict'
import { Error as ErrorComponent } from '@/components'

export type Metapath = ReturnType<FPL['toJSON']>

/**
 * Retreive output from the API, decode and return
 */
export function useMetapathOutput({ session_id, krg, head }: { session_id?: string, krg: KRG, head: Metapath }) {
  const { data: rawOutput, isLoading, error, mutate } = useSWRImmutable(() => head ? `${session_id ? `/api/socket/${session_id}` : ''}/api/db/process/${head.process.id}/output` : undefined)
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
  return { data: { output, outputNode }, error: error || decodeError, isLoading, mutate }
}

/**
 * Retreive inputs to this process from outputs of its inputs
 *  We rely on SWR to help us de-duplicate these requests
 */
export function useMetapathInputs({ session_id, krg, head }: { session_id?: string, krg: KRG, head: Metapath }) {
  const { data: rawInputs, error, isLoading } = useSWRMap<{ type: string, value: any }>(
    dict.values(head.process.inputs).map(({ id }) => `${session_id ? `/api/socket/${session_id}` : ''}/api/db/process/${id}/output`),
    fetcher
  )
  const processNode = krg.getProcessNode(head.process.type)
  const { inputs, decodeError } = React.useMemo(() => {
    try {
      if (!rawInputs) return { inputs: dict.isEmpty(processNode.inputs) ? {} as Record<string, unknown>: undefined }
      return {
        inputs: dict.init(
          dict.items(processNode.inputs).map(({ key, value }) => {
            if (Array.isArray(value)) {
              return {
                key,
                value: dict.items(head.process.inputs)
                  .filter(({ key: k }) => k.toString().startsWith(`${key}:`))
                  .map(({ value: { id } }) => {
                    const output = rawInputs[`${session_id ? `/api/socket/${session_id}` : ''}/api/db/process/${id}/output`]
                    if (!output) throw new Error(`No output for ${id}`)
                    if (output.type === ErrorComponent.spec) {
                      throw new Error(ErrorComponent.codec.decode(output.value))
                    }
                    return value[0].codec.decode(output.value)
                  })
              }
            } else {
              const output = rawInputs[`${session_id ? `/api/socket/${session_id}` : ''}/api/db/process/${head.process.inputs[key].id}/output`]
              if (!output) throw new Error(`No output for ${head.process.inputs[key].id}`)
              if (output.type === ErrorComponent.spec) {
                throw new Error(ErrorComponent.codec.decode(output.value))
              }
              return {
                key,
                value: value.codec.decode(output.value)
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
