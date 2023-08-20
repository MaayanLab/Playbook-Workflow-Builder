import React from 'react'
import dynamic from 'next/dynamic'
import Head from 'next/head'
import { view_report_icon } from '@/icons'
import Link from 'next/link'
import * as dict from '@/utils/dict'
import krg from '@/app/krg'
import { DataMetaNode } from '@/spec/metanode'
import { useRouter } from 'next/router'

import classNames from 'classnames'
import { PublicUserPlaybooks } from '@/app/api/client'
import { useAPIQuery } from '@/core/api/client'

const Icon = dynamic(() => import('@/app/components/icon'))
const Layout = dynamic(() => import('@/app/fragments/playbook/layout'))
const UserIdentity = dynamic(() => import('@/app/fragments/graph/useridentity'))

export default function CommunityPlaybooks() {
  const router = useRouter()
  const [inputFilters, setInputFilters] = React.useState<Record<string, boolean>>({})
  const [outputFilters, setOutputFilters] = React.useState<Record<string, boolean>>({})
  const [details, setDetails] = React.useState<Record<string, boolean>>({})
  const [search, setSearch] = React.useState('')
  const { data: playbooks } = useAPIQuery(PublicUserPlaybooks, {
    search,
    inputs: dict.items(inputFilters)
      .reduce((agg, { key, value }) => value ? [...agg, key] : agg, [] as string[])
      .join(', '),
    outputs: dict.items(outputFilters)
      .reduce((agg, { key, value }) => value ? [...agg, key] : agg, [] as string[])
      .join(', '),
  }, {
    keepPreviousData: true,
  })
  const { allInputs, allOutputs } = React.useMemo(() => {
    const allInputs: Record<string, DataMetaNode> = {}
    const allOutputs: Record<string, DataMetaNode> = {}
    if (playbooks) {
      playbooks.forEach(playbook => {
        (playbook.inputs||'').split(', ').forEach(input => {
          allInputs[input] = krg.getDataNode(input) as DataMetaNode
        });
        (playbook.outputs||'').split(', ').forEach(output => {
          allOutputs[output] = krg.getDataNode(output) as DataMetaNode
        });
      })
    }
    return { allInputs, allOutputs }
  }, [playbooks])
  return (
    <Layout>
      <Head>
        <title>Community Playbooks</title>
      </Head>

      <main className="flex-grow container mx-auto p-4 flex flex-col gap-6">
        <div className="flex flex-col gap-2">
          <div className="flex flex-row gap-2">
            <div className="flex-grow flex flex-col gap-2">
              <div className="prose"><h2>Inputs</h2></div>
              <div className="flex flex-row flex-wrap gap-2">
                {dict.values(allInputs).map(input => (
                  <button
                    key={input.spec}
                    onClick={() => {
                      setInputFilters(({ [input.spec]: cur, ...filters }) => {
                        if (!cur) return {...filters, [input.spec]: !cur}
                        else return filters
                      })
                    }}
                  >
                    <Icon
                      container="circle"
                      container_color={inputFilters[input.spec] ? '#B3CFFF' : '#ddd'}
                      size={1.5}
                      icon={input.meta.icon}
                      title={input.meta.label}
                    />
                  </button>
                ))}
              </div>
            </div>
            <div className="flex-grow flex flex-col gap-2">
              <div className="prose"><h2>Outputs</h2></div>
              <div className="flex flex-row flex-wrap gap-2">
                {dict.values(allOutputs).map(output => (
                  <button
                    key={output.spec}
                    onClick={() => {
                      setOutputFilters(({ [output.spec]: cur, ...filters }) => {
                        if (!cur) return {...filters, [output.spec]: !cur}
                        else return filters
                      })
                    }}
                  >
                    <Icon
                      key={output.spec}
                      container="circle"
                      container_color={outputFilters[output.spec] ? '#B3CFFF' : '#ddd'}
                      size={1.5}
                      icon={output.meta.icon}
                      title={output.meta.label}
                    />
                  </button>
                ))}
              </div>
            </div>
          </div>
        </div>
        <div className="bp4-input-group">
          <span className="bp4-icon bp4-icon-search" />
          <input
            type="search"
            className="bp4-input"
            placeholder="Search playbooks by title, description, and more"
            value={search}
            onChange={evt => {
              setSearch(() => evt.target.value)
            }}
          />
        </div>
        <div className="grid grid-rows-2 grid-cols-[min-content_max-content] md:grid-cols-6 md:grid-rows-[min-content_max-content]">
          <div className="bg-secondary font-bold p-3 text-center md:block hidden rounded-tl-lg"></div>
          <div className="bg-secondary font-bold p-3 text-center md:block hidden">Playbook</div>
          <div className="bg-secondary font-bold p-3 text-center md:block hidden">Inputs</div>
          <div className="bg-secondary font-bold p-3 text-center md:block hidden">Outputs</div>
          <div className="bg-secondary font-bold p-3 text-center md:block hidden">Clicks</div>
          <div className="bg-secondary font-bold p-3 text-center md:block hidden rounded-tr-lg">Actions</div>
          {playbooks ?
            playbooks.length > 0 ?
              playbooks.map(playbook =>
              <React.Fragment key={playbook.id}>
                <div className="col-span-1 row-span-1 p-3 hidden md:block">
                  <div className="tooltip z-50" data-tip="Click to view details">
                    <button
                      className="btn btn-ghost text-4xl"
                      onClick={evt => {
                        setDetails(({ [playbook.id]: cur, ...details }) => cur ? details : ({ ...details, [playbook.id]: true }))
                      }}
                    >{details[playbook.id] ? '-' : '+'}</button>
                  </div>
                </div>
                <div className={classNames('bg-secondary font-bold p-3 text-center md:hidden rounded-t-lg')}>Playbook</div>
                <div className="col-span-1 row-span-1 p-3">
                  <div className="prose md:h-12">{playbook.title}</div>
                </div>
                <div className={classNames('bg-secondary font-bold p-3 text-center md:hidden')}>Inputs</div>
                <div className="col-span-1 row-span-1 p-3">
                  <div className="flex flex-row flex-wrap gap-2 justify-center">
                    {(playbook.inputs||'').split(', ').map(spec => krg.getDataNode(spec)).map(input => (
                      <button
                        key={input.spec}
                        onClick={() => {
                          setInputFilters(({ [input.spec]: cur, ...filters }) => {
                            if (!cur) return {...filters, [input.spec]: !cur}
                            else return filters
                          })
                        }}
                      >
                        <Icon
                          container="circle"
                          container_color={inputFilters[input.spec] ? '#B3CFFF' : '#ddd'}
                          size={1.5}
                          icon={input.meta.icon}
                          title={input.meta.label}
                        />
                      </button>
                    ))}
                  </div>
                </div>
                <div className={classNames('bg-secondary font-bold p-3 text-center md:hidden')}>Outputs</div>
                <div className="col-span-1 row-span-1 p-3">
                  <div className="flex flex-row flex-wrap gap-2 justify-center">
                    {(playbook.outputs||'').split(', ').map(spec => krg.getDataNode(spec)).map(output => (
                      <button
                        key={output.spec}
                        onClick={() => {
                          setOutputFilters(({ [output.spec]: cur, ...filters }) => {
                            if (!cur) return {...filters, [output.spec]: !cur}
                            else return filters
                          })
                        }}
                      >
                        <Icon
                          key={output.spec}
                          container="circle"
                          container_color={outputFilters[output.spec] ? '#B3CFFF' : '#ddd'}
                          size={1.5}
                          icon={output.meta.icon}
                          title={output.meta.label}
                        />
                      </button>
                    ))}
                  </div>
                </div>
                <div className={classNames('bg-secondary font-bold p-3 text-center md:hidden')}>Clicks</div>
                <div className="col-span-1 row-span-1 p-3">
                  <div className="prose text-center md:h-12">
                    {playbook.clicks}
                  </div>
                </div>
                <div className={classNames('bg-secondary font-bold p-3 text-center md:hidden')}>Actions</div>
                <div className="col-span-1 row-span-1 p-3">
                  <div className="flex flex-row gap-2 justify-center">
                    <button onClick={() => {
                      // TODO: register click
                      router.push(`/report/${playbook.playbook}`)
                    }}>
                      <Icon icon={view_report_icon} color="black" title="Launch Playbook" />
                    </button>
                  </div>
                </div>
              {/* </div> */}
              <div className={classNames('bg-secondary font-bold p-3 row-span-1 col-span-1 text-center md:hidden rounded-b-lg')}>
                <button
                  className="btn btn-ghost text-4xl"
                  onClick={evt => {
                    setDetails(({ [playbook.id]: cur, ...details }) => cur ? details : ({ ...details, [playbook.id]: true }))
                  }}
                >{details[playbook.id] ? '-' : '+'}</button>
              </div>
              <div className={classNames('col-span-1 row-span-1 md:row-span-1 md:col-span-6')}>
                <div className={classNames('mx-auto prose whitespace-normal px-4', { hidden: !details[playbook.id] })}>
                  <p><b>Published</b>: {playbook.created.toString()}</p>
                  <p><b>Authors</b>:<br /><UserIdentity user={playbook.user} /></p>
                  <p><b>Description</b>: {playbook.description}</p>
                  <Link href={`/report/${playbook.playbook}`}><button className="bp4-button bp4-large">Launch</button></Link>
                </div>
              </div>
              <div className="col-span-2 row-span-1 md:hidden my-2">&nbsp;</div>
            </React.Fragment>
          )
          : <div className="col-span-2 row-span-1 md:row-span-1 md:col-span-6">
            <div className="alert">No playbooks currently registered matching this criteria.</div>
          </div>
          : null}
        </div>
        <Link href="/playbooks">
          <button className="btn btn-primary w-min whitespace-nowrap place-self-center text-center">
            Click for Playbooks Curated by the Playbook Partnership
          </button>
        </Link>
      </main>
    </Layout>
  )
}
