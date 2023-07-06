import React from 'react'
import dynamic from 'next/dynamic'
import Head from 'next/head'
import { view_report_icon } from '@/icons'
import Link from 'next/link'
import * as dict from '@/utils/dict'
import * as array from '@/utils/array'
import tsvector, { type TSVector } from "@/utils/tsvector"
import { DataMetaNode } from '@/spec/metanode'
import useAsyncEffect from 'use-async-effect'
import { useRouter } from 'next/router'

import ARCHS4_icon from '@/app/public/logos/datasources/ARCHS4.png'
import ChEA3_icon from '@/app/public/logos/datasources/ChEA3.png'
import ENCODE_icon from '@/app/public/logos/datasources/ENCODE.png'
import GTEx_icon from '@/app/public/logos/datasources/GTEx.png'
import GWAS_icon from '@/app/public/logos/datasources/GWAS.jpeg'
import IMPC_icon from '@/app/public/logos/datasources/IMPC.png'
import KOMP_icon from '@/app/public/logos/datasources/KOMP.png'
import LINCS_icon from '@/app/public/logos/datasources/LINCS.gif'
import IDG_icon from '@/app/public/logos/datasources/IDG.png'
import Metabolomics_icon from '@/app/public/logos/Metabolomics.png'
import KidsFirst_icon from '@/app/public/logos/KidsFirst.png'
import GlyGen_icon from '@/app/public/logos/glygen.svg'
import exRNA_icon from '@/app/public/logos/exRNA.png'
import Image, { StaticImageData } from 'next/image'
import classNames from 'classnames'

const dataSourceIcons: Record<string, StaticImageData> = {
  ARCHS4: ARCHS4_icon,
  ChEA: ChEA3_icon,
  ENCODE: ENCODE_icon,
  GTEx: GTEx_icon,
  'GWAS Catalog': GWAS_icon,
  IMPC: IMPC_icon,
  KOMP: KOMP_icon,
  'LINCS L1000': LINCS_icon,
  'Metabolomics': Metabolomics_icon,
  'GlyGen': GlyGen_icon,
  'KidsFirst': KidsFirst_icon,
  'exRNA': exRNA_icon,
  'IDG': IDG_icon,
}

const Icon = dynamic(() => import('@/app/components/icon'))
const Layout = dynamic(() => import('@/app/fragments/playbook/layout'))

type Playbook = {
  id: string,
  label: string,
  description: string,
  published: string,
  version: string,
  authors: string[],
  licenseUrl: string,
  license: string,
  url: string,
  dataSources: string[],
  inputs: Array<DataMetaNode>,
  outputs: Array<DataMetaNode>,
  clicks: number,
}

export default function Playbooks() {
  const router = useRouter()
  const [dataSourceFilters, setDataSourceFilters] = React.useState<Record<string, boolean>>({})
  const [inputFilters, setInputFilters] = React.useState<Record<string, boolean>>({})
  const [outputFilters, setOutputFilters] = React.useState<Record<string, boolean>>({})
  const [details, setDetails] = React.useState<Record<string, boolean>>({})
  const [playbooks, setPlaybooks] = React.useState<Array<Playbook>>()
  const [search, setSearch] = React.useState('')
  const playbook_tsvectors = React.useMemo(() => {
    const playbook_tsvectors: Record<string, TSVector> = {}
    playbooks?.forEach(playbook => {
      playbook_tsvectors[playbook.id] = tsvector([
        playbook.label,
        playbook.description,
        playbook.published,
        ...playbook.inputs.flatMap(input => [
          input.meta.label,
          input.meta.description,
        ]),
        ...playbook.outputs.flatMap(output => [
          output.meta.label,
          output.meta.description,
        ]),
        ...playbook.dataSources,
      ].join(' '))
    })
    return playbook_tsvectors
  }, [playbooks])
  const { allInputs, allOutputs, allDataSources, isLoading } = React.useMemo(() => {
    const allInputs: Record<string, DataMetaNode> = {}
    const allOutputs: Record<string, DataMetaNode> = {}
    const allDataSources: Record<string, string> = {}
    let isLoading = true
    if (playbooks) {
      playbooks.forEach(playbook => {
        playbook.inputs.forEach(input => {
          allInputs[input.spec] = input as DataMetaNode
        })
        playbook.outputs.forEach(output => {
          allOutputs[output.spec] = output as DataMetaNode
        })
        playbook.dataSources.forEach(dataSource => {
          allDataSources[dataSource] = dataSource
        })
      })
      isLoading = false
    }
    return { allInputs, allOutputs, allDataSources, isLoading }
  }, [playbooks])
  const filteredPlaybooks = React.useMemo(() => {
    if (!playbooks) return
    return playbooks.filter(playbook => {
      return array.all([
        ...dict.keys(inputFilters).map(filter => playbook.inputs.filter(input => input.spec === filter).length > 0),
        ...dict.keys(outputFilters).map(filter => playbook.outputs.filter(output => output.spec === filter).length > 0),
        ...dict.keys(dataSourceFilters).map(filter => playbook.dataSources.filter(dataSource => dataSource === filter).length > 0),
      ])
    })
  }, [playbooks, inputFilters, outputFilters, dataSourceFilters])
  const searchFilteredPlaybooks = React.useMemo(() => {
    if (!filteredPlaybooks) return
    if (!search) return filteredPlaybooks
    const search_tsvector = tsvector(search)
    const search_scores: Record<string, number> = {}
    filteredPlaybooks.forEach(playbook => {
      search_scores[playbook.id] = playbook_tsvectors[playbook.id]?.intersect(search_tsvector).size
    })
    const searchFilteredPlaybooks = filteredPlaybooks.filter(playbook => search_scores[playbook.id] > 0)
    searchFilteredPlaybooks.sort((a, b) => search_scores[b.id] - search_scores[a.id])
    return searchFilteredPlaybooks
  }, [filteredPlaybooks, search])
  useAsyncEffect(async (isMounted) => {
    const { default: demoPlaybooks } = await import('@/app/public/playbooksDemo')
    if (!isMounted()) return
    setPlaybooks(demoPlaybooks as Array<Playbook>)
  }, [])
  const DataSourceButton = React.useCallback(({ dataSource, size }: { dataSource: string, size: number }) => (
    <button
      key={dataSource}
      onClick={() => {
        setDataSourceFilters(({ [dataSource]: cur, ...filters }) => {
          if (!cur) return {...filters, [dataSource]: !cur}
          else return filters
        })
      }}
      className="flex flex-col place-items-center underline"
    >
      <span className={classNames('prose prose-sm', { 'text-shadow': dataSourceFilters[dataSource] })}>{dataSource}</span>
      {dataSource in dataSourceIcons ?
        <Image src={dataSourceIcons[dataSource]} objectFit="scale-down" width={size} height={size} />
        : null}
    </button>
  ), [dataSourceFilters, setDataSourceFilters])
  return (
    <Layout>
      <Head>
        <title>Playbooks</title>
      </Head>

      <main className="flex-grow container mx-auto p-4 flex flex-col gap-6">
        <div className="flex flex-col gap-2">
          <div className="flex flex-col gap-2">
            <div className="prose"><h2>Data Sources</h2></div>
            <div className="flex flex-row flex-wrap gap-2">
              {dict.values(allDataSources).map(dataSource => (
                <DataSourceButton
                  key={dataSource}
                  dataSource={dataSource}
                  size={80}
                />
              ))}
            </div>
          </div>
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
        <div className="grid grid-rows-2 grid-cols-[min-content_max-content] md:grid-cols-7 md:grid-rows-[min-content_max-content]">
          <div className="bg-secondary font-bold p-3 text-center md:block hidden rounded-tl-lg"></div>
          <div className="bg-secondary font-bold p-3 text-center md:block hidden">Playbook</div>
          <div className="bg-secondary font-bold p-3 text-center md:block hidden">Inputs</div>
          <div className="bg-secondary font-bold p-3 text-center md:block hidden">Outputs</div>
          <div className="bg-secondary font-bold p-3 text-center md:block hidden">Sources</div>
          <div className="bg-secondary font-bold p-3 text-center md:block hidden">Clicks</div>
          <div className="bg-secondary font-bold p-3 text-center md:block hidden rounded-tr-lg">Actions</div>
          {searchFilteredPlaybooks ?
            searchFilteredPlaybooks.length > 0 ?
              searchFilteredPlaybooks.map(playbook =>
              <React.Fragment key={playbook.id}>
                <div className="col-span-1 row-span-1 p-3 hidden md:block">
                  <div className="tooltip z-50" data-tip="Click to view details">
                    <button
                      className="btn btn-ghost text-4xl text-black dark:text-white"
                      onClick={evt => {
                        setDetails(({ [playbook.id]: cur, ...details }) => cur ? details : ({ ...details, [playbook.id]: true }))
                      }}
                    >{details[playbook.id] ? '-' : '+'}</button>
                  </div>
                </div>
                <div className={classNames('bg-secondary font-bold p-3 text-center md:hidden rounded-t-lg')}>Playbook</div>
                <div className="col-span-1 row-span-1 p-3">
                  <div className="prose md:h-12">{playbook.label}</div>
                </div>
                <div className={classNames('bg-secondary font-bold p-3 text-center md:hidden')}>Inputs</div>
                <div className="col-span-1 row-span-1 p-3">
                  <div className="flex flex-row flex-wrap gap-2 justify-center">
                    {playbook.inputs.map(input => (
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
                    {playbook.outputs.map(output => (
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
                <div className={classNames('bg-secondary font-bold p-3 text-center md:hidden')}>Sources</div>
                <div className="col-span-1 row-span-1 p-3">
                  <div className="flex flex-row flex-wrap gap-2 justify-center">
                    {playbook.dataSources.map(dataSource => (
                      <DataSourceButton
                        key={dataSource}
                        dataSource={dataSource}
                        size={50}
                      />
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
                      router.push(playbook.url)
                    }}>
                      <Icon icon={view_report_icon} className="fill-black dark:fill-white" title="Launch Playbook" />
                    </button>
                  </div>
                </div>
              {/* </div> */}
              <div className={classNames('bg-secondary font-bold p-3 row-span-1 col-span-1 text-center md:hidden rounded-b-lg')}>
                <button
                  className="btn btn-ghost text-4xl text-black dark:text-white"
                  onClick={evt => {
                    setDetails(({ [playbook.id]: cur, ...details }) => cur ? details : ({ ...details, [playbook.id]: true }))
                  }}
                >{details[playbook.id] ? '-' : '+'}</button>
              </div>
              <div className={classNames('col-span-1 row-span-1 md:row-span-1 md:col-span-7')}>
                <div className={classNames('mx-auto prose whitespace-normal px-4', { hidden: !details[playbook.id] })}>
                  <div className="flex flex-row gap-2 flex-wrap">
                    <div className="badge badge-primary">v{playbook.version}</div>
                    <div className="badge bg-gray-300 dark:bg-gray-600 border-none"><a className="no-underline text-blue-700" href={playbook.licenseUrl} target="blank">{playbook.license}</a></div>
                  </div>
                  <p><b>Published</b>: {playbook.published}</p>
                  <p><b>Authors</b>:<br />{playbook.authors.join(', ')}</p>
                  <p><b>Description</b>: {playbook.description}</p>
                  <Link href={playbook.url}><button className="bp4-button bp4-large">Launch</button></Link>
                </div>
              </div>
              <div className="col-span-2 row-span-1 md:hidden my-2">&nbsp;</div>
            </React.Fragment>
          )
          : <div className="col-span-7">
            <div className="alert prose">No playbooks currently registered matching this criteria.</div>
          </div>
          : null}
        </div>
        <Link href="/playbooks/community">
          <button className="btn btn-primary w-min whitespace-nowrap place-self-center text-center">
            Click for Additional Playbooks Contributed by Other Users
          </button>
        </Link>
      </main>
    </Layout>
  )
}
