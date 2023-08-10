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
import { clientLoadExample } from  '@/components/core/file/api/example.h5ad/client'

const Bp4Button = dynamic(() => import('@blueprintjs/core').then(({ Button }) => Button))

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
      <h2>File: {url}</h2>
    </div>
  ))
  .build()

export const FileInput = MetaNode('FileInput')
  .meta({
    label: 'Input File',
    description: 'Upload a Data File',
    icon: [input_icon],
  })
  .inputs()
  .output(FileURL)
  .prompt(props => {
    const [currentFile, setCurrentFile] = React.useState<{ description?: string | null, url?: string, filename?: string, size?: number }>({})
    const [tab, setTab] = React.useState<'upload' | 'url'>()
    const fileInputRef = React.useRef<HTMLInputElement>(null)
    const { data: session } = useSessionWithId()
    const output = React.useMemo(() => props.output || { description: undefined, url: undefined, filename: undefined, size: undefined }, [props.output])
    React.useEffect(() => {
      setCurrentFile(output)
    }, [output])
    const activeTab = tab ? tab : currentFile.url ? 'url' : 'upload'
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
              if (currentFile.url && currentFile.url === props.output?.url && currentFile.filename && currentFile.size) {
                props.submit({
                  description: currentFile.description,
                  url: currentFile.url,
                  filename: currentFile.filename,
                  size: currentFile.size,
                })
              } else if (activeTab === 'upload') {
                const formData = new FormData(evt.currentTarget)
                const rawDescription = formData.get('description')
                const description = rawDescription === null ? undefined : rawDescription.toString()
                const { file: [record] } = await clientUploadFile(formData)
                props.submit({ description, ...record })
              } else if (currentFile.url) {
                const formData = new FormData(evt.currentTarget)
                const rawDescription = formData.get('description')
                const description = rawDescription === null ? undefined : rawDescription.toString()
                const record = await clientFetchFile({ url: currentFile.url })
                props.submit({ description, ...record })
              } else {
                console.error('Missing input')
              }
            }}
          >
            <div className="tabs">
              <button
                type="button"
                className={classNames("tab tab-bordered", { "tab-active": activeTab === 'upload' })}
                onClick={() => { setTab('upload') }}
              >Upload file</button>
              <button
                type="button"
                className={classNames("tab tab-bordered", { "tab-active": activeTab === 'url' })}
                onClick={() => { setTab('url') }}
              >File from URL</button>
            </div>
            <div className="h-24 flex flex-col items-stretch justify-center">
              <div className={classNames({ 'hidden': activeTab !== 'upload' })}>
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
              <div className={classNames({ 'hidden': activeTab !== 'url' })}>
                <div className="bp4-input-group">
                  <input
                    type="text"
                    className="bp4-input w-full"
                    placeholder="e.g. drs://somehost/id"
                    value={currentFile.url}
                    onChange={evt => {
                      setCurrentFile(({ url: _, ...file }) => ({ ...file, url: evt.target.value }))
                    }}
                  />
                </div>
              </div>
            </div>
            <div className="bp4-input-group">
              <input
                type="text"
                name="description"
                className="bp4-input"
                placeholder={`File description`}
                onChange={evt => {setCurrentFile(({ description: _, ...file }) => ({ ...file, description: evt.target.value }))}}
                value={currentFile.description||''}
              />
            </div>
            <div className="inline-flex flex-row">
              <Bp4Button
                large
                text="Example"
                rightIcon="bring-data"
                onClick={async () => {
                  props.submit(await clientLoadExample())
                }}
              />
              <Bp4Button
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
  })
  .story(props =>
    `A file${props.output?.description ? ` containing ${props.output?.description}` : ''} was first uploaded.`
  )
  .build()
