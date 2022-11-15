import React from 'react'
import type { GetServerSidePropsContext } from 'next'
import { useRouter } from 'next/router'
import fpprg from '@/app/fpprg'
import { FPL } from '@/core/FPPRG'
import krg from '@/app/krg'
import { z } from 'zod'
import useSWRImmutable from 'swr/immutable'
import { SWRConfig } from 'swr'
import dynamic from 'next/dynamic'
import type KRG from '@/core/KRG'
import Head from 'next/head'

const Header = dynamic(() => import('@/app/fragments/playbook/header'))
const Footer = dynamic(() => import('@/app/fragments/playbook/footer'))
const Button = dynamic(() => import('@blueprintjs/core').then(({ Button }) => Button))

type Metapath = ReturnType<FPL['toJSON']>

const QueryType = z.object({
  id: z.string().optional(),
})

export async function getServerSideProps(ctx: GetServerSidePropsContext) {
  const params = QueryType.parse(ctx.params || {})
  if (!params.id) {
    return {
      props: {
        fallback: {}
      }
    }
  }
  const fpl = fpprg.getFPL(params.id)
  if (fpl === undefined) {
    return { notFound: true }
  }
  const results = fpl.resolve().map(fpl => fpl.toJSON())
  const fallback: Record<string, unknown> = {
    [`/api/db/fpl/${params.id}`]: results
  }
  for (const result of results) {
    const output = fpprg.getResolved(result.process.id)
    if (output) fallback[`/api/db/process/${result.process.id}/output`] = output.toJSON().data
  }
  return {
    props: {
      fallback,
    }
  }
}

async function fetcher(path: string): Promise<Array<Metapath>> {
  const req = await fetch(path)
  return await req.json()
}

function useSticky(value: any) {
  const ref = React.useRef()
  if (ref !== undefined && value !== undefined) ref.current = value
  return ref.current
}

function useSWRImmutableSticky<T = any>(key?: string) {
  const swr = useSWRImmutable<T>(key)
  const stickyData = useSticky(swr.data)
  return { ...swr, data: swr.data === undefined ? stickyData : swr.data }
}

export default function App({ fallback }: { fallback: any }) {
  const router = useRouter()
  const params = QueryType.parse(router.query)
  return (
    <div className="flex flex-col min-w-screen min-h-screen">
      <Head>
        <title>Playbook</title>
      </Head>

      <Header />

      <SWRConfig value={{ fallback, fetcher }}>
        <main className="flex-grow container mx-auto py-4 flex flex-col">
          <Cells krg={krg} id={params.id} />
        </main>
      </SWRConfig>

      <Footer />
    </div>
  )
}

function Cells({ krg, id }: { krg: KRG, id?: string }) {
  const router = useRouter()
  const { data: metapath, error } = useSWRImmutableSticky<Array<Metapath>>(id ? `/api/db/fpl/${id}` : undefined)
  const head = metapath ? metapath[metapath.length - 1] : undefined
  const processNode = head ? krg.getProcessNode(head.process.type) : undefined
  const actions = processNode ? krg.getNextProcess(processNode.output.spec) : krg.getNextProcess('')
  return (
    <div className="container mx-auto py-4">
      {error ? <div>{error}</div> : null}
      {(metapath||[]).map((head, index) =>
        <Cell key={index} krg={krg} id={id} head={head} />
      )}
      {actions.length > 0 ? (
        <div className="flex-grow flex-shrink bp4-card my-2">
          <h2 className="bp4-heading">Actions</h2>
          {actions.map(proc =>
            <div key={proc.spec}>
              {Object.keys(proc.inputs).length > 0 ? (
                <>
                  <span className="bg-secondary rounded-full p-3">{Object.values(proc.inputs).map((i) => i.spec).join(', ')}</span>
                  <span> =&gt; </span>
                </>
              ) : null}
              <Button
                large
                onClick={async () => {
                  const inputs: Record<string, { id: string }> = {}
                  if (head) {
                    for (const i in proc.inputs) {
                      inputs[i] = { id: head.process.id }
                    }
                  }
                  const req = await fetch(`/api/db/fpl/${id || 'start'}/extend`, {
                    method: 'POST',
                    body: JSON.stringify({
                      type: proc.spec,
                      inputs,
                    })
                  })
                  const res = z.string().parse(await req.json())
                  router.push(`/report/${res}`, undefined, { shallow: true, scroll: false })
                }}
              >{proc.spec}</Button>
              <span> =&gt; </span>
              <span className="bg-secondary rounded-full p-3">{proc.output.spec}</span>
            </div>
          )}
        </div>
      ) : null}
    </div>
  )
}


function Cell({ krg, id, head }: { krg: KRG, id?: string, head: Metapath }) {
  const router = useRouter()
  const { data: rawOutput, error: outputError } = useSWRImmutableSticky<any>(`/api/db/process/${head.process.id}/output`)
  const processNode = krg.getProcessNode(head.process.type)
  const inputs: any = head.process.inputs
  const outputNode = rawOutput && !outputError ? krg.getDataNode(rawOutput.type) : processNode.output
  const output = rawOutput && outputNode ? outputNode.codec.decode(rawOutput.value) : rawOutput
  const View = outputNode ? outputNode.view : undefined
  const Prompt = 'prompt' in processNode ? processNode.prompt : undefined
  return (
    <div className="flex-grow flex flex-col">
      <div className="flex-grow flex-shrink bp4-card my-2">
        <h2 className="bp4-heading">{processNode.meta.label || processNode.spec}</h2>
        {Prompt ? <Prompt
          inputs={inputs}
          output={output}
          submit={async (output) => {
            const req = await fetch(`/api/db/fpl/${id}/rebase/${head.process.id}`, {
              method: 'POST',
              body: JSON.stringify({
                type: head.process.type,
                data: {
                  type: processNode.output.spec,
                  value: processNode.output.codec.encode(output),
                },
                inputs,
              })
            })
            const res = z.object({ head: z.string(), rebased: z.string() }).parse(await req.json())
            router.push(`/report/${res.head}`, undefined, { shallow: true, scroll: false })
          }}
        />
        : processNode.meta.description ? <p className="bp4-ui-text">{processNode.meta.description}</p>
        : null}
      </div>
      <div className="flex-grow py-4 bp4-card">
        {outputNode ? (
          <>
            <h2 className="bp4-heading">{outputNode.meta.label || outputNode.spec}</h2>
            {View && output ? View(output) : 'Waiting for input'}
          </>
        ) : (
          <div>Loading...</div>
        )}
      </div>
    </div>
  )
}
