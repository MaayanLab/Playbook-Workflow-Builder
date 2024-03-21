import React from 'react'
import { MetaNode } from '@/spec/metanode'
import { z } from 'zod'
import dynamic from 'next/dynamic'
import { file_icon, input_icon } from '@/icons'
import * as Auth from 'next-auth/react'
import { useSessionWithId } from '@/app/extensions/next-auth/hooks'
import classNames from 'classnames'
import { clientUploadFile } from  '@/components/core/file/api/upload/client'
import { clientFetchFile } from  '@/components/core/file/api/fetch/client'

const Bp5Button = dynamic(() => import('@blueprintjs/core').then(({ Button }) => Button))

export const FileC = z.object({
  description: z.string().optional().nullable(),
  url: z.string(),
  filename: z.string(),
  size: z.number().optional(),
})

export const FileURL = MetaNode('FileURL')
  .meta({
    label: 'File URL',
    description: 'URL to a File',
    icon: [file_icon],
  })
  .codec(FileC)
  .view(({ url }) => (
    <div>
      <h2 className="prose">File: {url}</h2>
    </div>
  ))
  .build()

export function FilePrompt(props: {
  example?(): Promise<{
    url: string;
    filename: string;
    sha256: string;
    size: number;
  }>,
  data?: {
      url: string;
      filename: string;
      description?: string | null | undefined;
      size?: number | undefined;
  } | undefined;
  inputs: {};
  submit: (data: {
      url: string;
      filename: string;
      description?: string | null | undefined;
      size?: number | undefined;
  }) => void;
  session_id?: string | undefined;
}) {
  const [currentFile, setCurrentFile] = React.useState<{ description?: string | null, url?: string, filename?: string, size?: number }>({})
  const [tab, setTab] = React.useState<'upload' | 'url'>('upload')
  const fileInputRef = React.useRef<HTMLInputElement>(null)
  const { data: session } = useSessionWithId()
  const output = React.useMemo(() => props.data || { description: undefined, url: undefined, filename: undefined, size: undefined }, [props.data])
  React.useEffect(() => {
    setCurrentFile(output)
    if (output && output.url && !output.url.startsWith((process.env.NEXT_PUBLIC_URL as string).replace(/https?:/, 'drs:'))) {
      // only go to url when the url isn't our DRS server
      setTab('url')
    }
  }, [output])
  return (
    <div>
      {!session || !session.user ? (
        <div className="alert alert-warning shadow-lg block prose">
          You are required to &nbsp; <button className="btn btn-sm" onClick={() => {Auth.signIn()}}>sign in</button> &nbsp; to upload files.
        </div>
      ) : (
        <form
          className="my-2 inline-flex flex-col w-full max-w-sm"
          onSubmit={async (evt) => {
            evt.preventDefault()
            if (currentFile.url && currentFile.url === props.data?.url && currentFile.filename && currentFile.size) {
              props.submit({
                description: currentFile.description,
                url: currentFile.url,
                filename: currentFile.filename,
                size: currentFile.size,
              })
            } else if (tab === 'upload') {
              const formData = new FormData(evt.currentTarget)
              const rawDescription = formData.get('description')
              const description = rawDescription === null ? undefined : rawDescription.toString()
              const { file: [record] } = await clientUploadFile(formData, props.session_id)
              props.submit({ description, ...record })
            } else if (currentFile.url) {
              const formData = new FormData(evt.currentTarget)
              const rawDescription = formData.get('description')
              const description = rawDescription === null ? undefined : rawDescription.toString()
              const record = await clientFetchFile({ url: currentFile.url }, props.session_id)
              props.submit({ description, ...record })
            } else {
              console.error('Missing input')
            }
          }}
        >
          <div className="tabs">
            <button
              type="button"
              className={classNames("tab tab-bordered", { "tab-active": tab === 'upload' })}
              onClick={() => { setTab('upload') }}
            >Upload file</button>
            <button
              type="button"
              className={classNames("tab tab-bordered", { "tab-active": tab === 'url' })}
              onClick={() => { setTab('url') }}
            >File from URL</button>
          </div>
          <div className="h-24 flex flex-col items-stretch justify-center">
            <div className={classNames({ 'hidden': tab !== 'upload' })}>
              <input
                ref={fileInputRef}
                name="file"
                className={classNames('file-input file-input-lg', { 'hidden': currentFile.url })}
                onChange={evt => {
                  const url = evt.currentTarget && evt.currentTarget.files && evt.currentTarget.files[0] ? evt.currentTarget.files[0].name : output.url
                  setCurrentFile(({ description, ...file }) => ({ description, url }))
                }}
                type="file"
              />
              <div
                className={classNames('inline-flex flex-row items-center gap-4', { 'hidden': !currentFile.url })}
                onClick={evt => {
                  evt.preventDefault()
                  if (fileInputRef.current) {
                    fileInputRef.current.click()
                  }
                }}
              >
                <button className="btn btn-lg">Choose File</button>
                <span className="text-lg">{currentFile.filename||currentFile.url||'No file chosen'}</span>
              </div>
            </div>
            <div className={classNames({ 'hidden': tab !== 'url' })}>
              <div className="bp5-input-group">
                <input
                  type="text"
                  className="bp5-input w-full"
                  placeholder="e.g. drs://somehost/id"
                  value={currentFile.url}
                  onChange={evt => {
                    setCurrentFile(({ url: _, ...file }) => ({ ...file, url: evt.target.value }))
                  }}
                />
              </div>
            </div>
          </div>
          <div className="bp5-input-group">
            <input
              type="text"
              name="description"
              className="bp5-input"
              placeholder={`File description`}
              onChange={evt => {setCurrentFile(({ description: _, ...file }) => ({ ...file, description: evt.target.value }))}}
              value={currentFile.description||''}
            />
          </div>
          <div className="inline-flex flex-row">
            {props.example ?
              <Bp5Button
                large
                text="Example"
                rightIcon="bring-data"
                onClick={async () => {
                  if (!props.example) return
                  props.submit(await props.example())
                }}
              />
              : null}
            <Bp5Button
              large
              disabled={
                currentFile.url === undefined
                || (
                  currentFile.description === output.description
                  && currentFile.url === output.url
                )
              }
              type="submit"
              text="Submit"
              rightIcon="send-to-graph"
            />
          </div>
        </form>
      )}
    </div>
  )
}

export const FileInput = MetaNode('FileInput')
  .meta({
    label: 'Input File',
    description: 'Upload a Data File',
    tags: {
      Type: {
        File: 1,
      },
      Cardinality: {
        Matrix: 1,
      },
    },
    icon: [input_icon],
    hidden: true,
  })
  .inputs()
  .output(FileURL)
  .prompt(FilePrompt)
  .story(props =>
    `A file${props.output?.description ? ` containing ${props.output?.description}` : ''} was first uploaded.`
  )
  .build()
