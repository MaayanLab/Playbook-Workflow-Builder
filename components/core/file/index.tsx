import React from 'react'
import { MetaNode } from '@/spec/metanode'
import { z } from 'zod'
import dynamic from 'next/dynamic'
import { file_icon, input_icon } from '@/icons'
import * as Auth from 'next-auth/react'
import { useSessionWithId } from '@/app/extensions/next-auth/hooks'

const Bp4Button = dynamic(() => import('@blueprintjs/core').then(({ Button }) => Button))

export const FileURL = MetaNode('FileURL')
  .meta({
    label: 'File URL',
    description: 'URL to a File',
    icon: [file_icon],
  })
  .codec(z.object({ description: z.string().optional(), url: z.string() }))
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
    default: '',
    icon: [input_icon],
  })
  .inputs()
  .output(FileURL)
  .prompt(props => {
    const [currentFile, setCurrentFile] = React.useState<{ description?: string, url?: string }>({})
    const fileInputRef = React.useRef<HTMLInputElement>(null)
    const { data: session } = useSessionWithId()
    const output = React.useMemo(() => props.output || { description: undefined, url: undefined }, [props.output])
    React.useEffect(() => {
      setCurrentFile(output)
    }, [output])
    return (
      <div>
        {!session || !session.user ? (
          <div className="alert alert-warning shadow-lg block">
            You are required to &nbsp; <button className="btn btn-sm" onClick={() => {Auth.signIn()}}>sign in</button> &nbsp; to upload files.
          </div>
        ) : (
          <form
            className="my-2 inline-flex flex-col"
            onSubmit={async (evt) => {
              evt.preventDefault()
              if (currentFile.url && currentFile.url === props.output?.url) {
                props.submit({ description: currentFile.description, url: currentFile.url })
              } else {
                const formData = new FormData(evt.currentTarget)
                const rawDescription = formData.get('description')
                const description = rawDescription === null ? undefined : rawDescription.toString()
                const res = await fetch(`/api/components/core/file/upload`, { method: 'POST', body: formData })
                const records: { file: string[] } = await res.json()
                props.submit({ description, url: records.file[0] })
              }
            }}
          >
            <input
              ref={fileInputRef}
              name="file"
              className={`file-input file-input-lg ${currentFile.url ? 'hidden' : ''}`}
              onChange={evt => {
                const url = evt.currentTarget && evt.currentTarget.files && evt.currentTarget.files[0] ? evt.currentTarget.files[0].name : output.url
                setCurrentFile(({ description, url: _ }) => ({ description, url }))
              }}
              type="file"
            />
            <div
              className={`inline-flex flex-row items-center gap-4 ${currentFile.url ? '' : 'hidden'}`}
              onClick={evt => {
                evt.preventDefault()
                if (fileInputRef.current) {
                  fileInputRef.current.click()
                }
              }}
            >
              <button className="btn btn-lg">Choose File</button>
              <span className="text-lg">{currentFile.url||'No file chosen'}</span>
            </div>
            <div className="bp4-input-group">
              <input
                type="text"
                name="description"
                className="bp4-input"
                placeholder={`File description`}
                onChange={evt => {setCurrentFile(({ description: _, url }) => ({ url, description: evt.target.value }))}}
                value={currentFile.description||''}
              />
            </div>
            <div className="inline-flex flex-row">
              <a
                href="/api/components/core/file/example_matrix.tsv"
                download="example_matrix.tsv"
              >
                <Bp4Button
                  large
                  text="Example"
                  rightIcon="bring-data"
                />
              </a>
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
