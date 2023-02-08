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
import Link from 'next/link'
import Icon from '@/app/components/icon'
import { status_awaiting_input_icon, status_complete_icon, status_waiting_icon, status_alert_icon, view_in_graph_icon, fork_icon, share_icon, start_icon, func_icon, variable_icon } from '@/icons'
import * as array from '@/utils/array'
import * as dict from '@/utils/dict'
import usePublicUrl from '@/utils/next-public-url'

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
  const fpl = await fpprg.getFPL(params.id)
  if (fpl === undefined) {
    return { notFound: true }
  }
  const results = fpl.resolve().map(fpl => fpl.toJSON())
  const fallback: Record<string, unknown> = {
    [`/api/db/fpl/${params.id}`]: results
  }
  for (const result of results) {
    const output = await fpprg.getResolved(result.process.id)
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
  return { ...swr, data: swr.data === undefined ? stickyData : swr.data, isLoading: swr.data === undefined }
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

function ShareButton({ id }: { id: string | undefined }) {
  const publicUrl = usePublicUrl({ absolute: true })
  const [share, setShare] = React.useState(false)
  const onClick = React.useCallback(() => {
    const graphUrl = document.getElementById('graph-url') as HTMLInputElement
    graphUrl.select()
    graphUrl.setSelectionRange(0, 99999)
    navigator.clipboard.writeText(graphUrl.value)
  }, [id])
  React.useEffect(() => { if (share) { onClick() } }, [share])
  return (
    <>
      <button className={`bp4-button bp4-minimal${share ? ' hidden' : ''}`} onClick={() => {setShare(true)}}>
        <Icon icon={share_icon} color="black" />
      </button>
      <div className={`bp4-control-group inline-block${share ? '': ' hidden'}`}>
        <input id="graph-url" type="text" className="bp4-input" value={`${publicUrl}/report${id ? `/${id}` : ''}`} readOnly />
        <button className="bp4-button bp4-icon-link" onClick={onClick} />
        <button className="bp4-button bp4-icon-cross" onClick={() => {setShare(false)}} />
      </div>
    </>
  )
}

function Cells({ krg, id }: { krg: KRG, id?: string }) {
  const router = useRouter()
  const { data: metapath, error } = useSWRImmutableSticky<Array<Metapath>>(id ? `/api/db/fpl/${id}` : undefined)
  const head = metapath ? metapath[metapath.length - 1] : undefined
  const selections = metapath ? dict.init(metapath.map(item => ({ key: item.process.id, value: { process: item.process, processNode: krg.getProcessNode(item.process.type) } }))) : {}
  const processNode = head ? krg.getProcessNode(head.process.type) : undefined
  const actions = processNode ? krg.getNextProcess(processNode.output.spec) : krg.getNextProcess('')
  return (
    <div className="flex flex-col py-4 gap-2">
      <div className="flex-grow flex-shrink bp4-card p-0">
        <div className="p-3">
          <div className="flex flex-row gap-2">
            <Icon icon={start_icon} />
            <h2 className="bp4-heading">
              Playbook
            </h2>
          </div>
        </div>
        {error ? <div className="alert alert-error">{error}</div> : null}
        <div className="border-t-secondary border-t-2 mt-2">
          <Link href={`/graph${id ? `/${id}/node/start` : ``}`}>
            <button className="bp4-button bp4-minimal">
              <Icon icon={view_in_graph_icon} />
            </button>
          </Link>
          <Link href={`/graph${id ? `/${id}/node/start/extend` : `/start/extend`}`}>
            <button className="bp4-button bp4-minimal">
              <Icon icon={fork_icon} color="black" />
            </button>
          </Link>
          <ShareButton id={id} />
        </div>
      </div>
      {(metapath||[]).map((head, index) =>
        <Cell key={index} krg={krg} index={index} id={id} head={head} />
      )}
      {actions.length > 0 ? (
        <div className="flex-grow flex-shrink bp4-card p-0">
          <div className="p-4">
            <h2 className="bp4-heading">Actions</h2>
            {actions.map(proc => {
              // determine if multi-inputs are satisfiable
              const disabled = !array.all(
                dict.values(proc.inputs).map((value) => {
                  if (Array.isArray(value)) {
                    return dict.values(selections).filter(selection => selection.processNode.output.spec === value[0].spec).length > 1
                  } else {
                    return dict.values(selections).filter(selection => selection.processNode.output.spec === value.spec).length >= 1
                  }
                })
              )
              if (disabled) return null
              return (
                <div key={proc.spec}>
                  {!dict.isEmpty(proc.inputs) ? (
                    <>
                      <span className="bg-secondary rounded-full p-3">{dict.values(proc.inputs).map((i) => Array.isArray(i) ? `[${array.ensureOne(i).meta.label}]` : i.meta.label).join(', ')}</span>
                      <span> =&gt; </span>
                    </>
                  ) : null}
                  <Button
                    large
                    onClick={async () => {
                      const inputs: Record<string, { id: string }> = {}
                      dict.items(proc.inputs).forEach(({ key: arg, value: input }) => {
                        if (Array.isArray(input)) {
                          dict.values(selections)
                            .filter(selection => selection.processNode.output.spec === input[0].spec)
                            .forEach((selection, i) => {
                              inputs[`${arg as string}:${i}`] = { id: selection.process.id }
                            })
                        } else {
                          dict.values(selections).forEach(selection => {
                            if (selection.processNode.output.spec === input.spec) {
                              inputs[arg as string] = { id: selection.process.id }
                            }
                          })
                        }
                      })
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
                  >{proc.meta.label}</Button>
                  <span> =&gt; </span>
                  <span className="bg-secondary rounded-full p-3">{proc.output.meta.label}</span>
                </div>
              )
            })}
          </div>
          <div className="border-t-secondary border-t-2 mt-2">
            <Link href={`/graph${id ? `/${id}/node/${id}/extend` : `/start/extend`}`}>
              <button className="bp4-button bp4-minimal">
                <Icon icon={view_in_graph_icon} />
              </button>
            </Link>
          </div>
        </div>
      ) : null}
    </div>
  )
}

function Cell({ krg, index, id, head }: { krg: KRG, index: number, id?: string, head: Metapath }) {
  const router = useRouter()
  const { data: rawOutput, error: outputError, isLoading } = useSWRImmutableSticky<any>(`/api/db/process/${head.process.id}/output`)
  const processNode = krg.getProcessNode(head.process.type)
  const inputs: any = head.process.inputs
  const outputNode = rawOutput && !outputError ? krg.getDataNode(rawOutput.type) : processNode.output
  const output = rawOutput && outputNode ? outputNode.codec.decode(rawOutput.value) : rawOutput
  const View = outputNode ? outputNode.view : undefined
  const Prompt = 'prompt' in processNode ? processNode.prompt : undefined
  return (
    <>
      <div className="flex-grow flex-shrink items-center overflow-auto bp4-card p-0">
        <div className="p-3">
          <div className="flex flex-row gap-2">
            <Icon icon={processNode.meta.icon || func_icon} />
            <h2 className="bp4-heading">{processNode.meta.label || processNode.spec}</h2>
          </div>
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
        <div className="border-t-secondary border-t-2 mt-2">
          <Link href={`/graph/${id}/node/${head.id}`}>
            <button className="bp4-button bp4-minimal">
              <Icon icon={view_in_graph_icon} />
            </button>
          </Link>
        </div>
      </div>
      <div className="flex-grow flex-shrink items-center overflow-auto bp4-card p-0">
        <div className="p-3">
          <div className="flex flex-row gap-2">
            <Icon icon={(outputNode && outputNode.meta.icon) || variable_icon} />
            <h2 className="bp4-heading">{(outputNode && (outputNode.meta.label || processNode.spec)) || "Loading"}</h2>
          </div>
          {outputNode && View && output ? View(output) : isLoading ? 'Waiting for results' : 'Waiting for input'}
        </div>
        <div className="border-t-secondary border-t-2 mt-2">
          <Link href={`/graph/${id}/node/${head.id}`}>
            <button className="bp4-button bp4-minimal">
              <Icon icon={view_in_graph_icon} />
            </button>
          </Link>
          <Link href={`/graph/${id}/node/${head.id}/extend`}>
            <button className="bp4-button bp4-minimal">
              <Icon icon={fork_icon} color="black" />
            </button>
          </Link>
          <button className="bp4-button bp4-minimal" disabled>
            {isLoading ?
              <Icon icon={status_waiting_icon} color="#DAA520" />
              : (outputNode ?
                  (output ?
                    (outputNode.spec === 'Error' ?
                      <Icon icon={status_alert_icon} color="#DC143C" />
                      : <Icon icon={status_complete_icon} color="#008000" />)
                    : <Icon icon={status_awaiting_input_icon} color="#B8860B" />)
                  : <Icon icon={status_waiting_icon} color="#DAA520" />)}
          </button>
        </div>
      </div>
    </>
  )
}
