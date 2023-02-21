import React from 'react'
import { MetaNode } from '@/spec/metanode'
import { z } from 'zod'
import dynamic from 'next/dynamic'
import { file_icon, input_icon } from '@/icons'
import * as Auth from 'next-auth/react'
import { useSessionWithId } from '@/app/extensions/next-auth/hooks'

const Bp4Button = dynamic(() => import('@blueprintjs/core').then(({ Button }) => Button))

export const FileURL = MetaNode.createData('FileURL')
  .meta({
    label: 'File URL',
    description: 'An arbitrary file url',
    icon: [file_icon],
  })
  .codec(z.string())
  .view(file => (
    <div>
      <h2>File: {file}</h2>
    </div>
  ))
  .build()

export const FileInput = MetaNode.createProcess('FileInput')
  .meta({
    label: 'Input a File',
    description: 'A file upload',
    default: '',
    icon: [input_icon],
  })
  .inputs()
  .output(FileURL)
  .prompt(props => {
    const [currentFile, setCurrentFile] = React.useState<boolean>(false)
    const fileInputRef = React.useRef<HTMLInputElement>(null)
    const { data: session } = useSessionWithId()
    return (
      <div>
        {session === null ? (
          <div className="alert alert-warning shadow-lg block">
            You are required to &nbsp; <button className="btn btn-sm" onClick={() => {Auth.signIn()}}>sign in</button> &nbsp; to upload files.
          </div>
        ) : (
          <form
            className="my-2 inline-flex flex-col"
            onSubmit={async (evt) => {
              evt.preventDefault()
              const formData = new FormData(evt.currentTarget)
              const res = await fetch(`/api/components/core/file/upload`, { method: 'POST', body: formData })
              const records: { file: string[] } = await res.json()
              props.submit(records.file[0])
            }}
          >
            <input
              ref={fileInputRef}
              name="file"
              className={`file-input file-input-lg ${props.output && !currentFile ? 'hidden' : ''}`}
              onChange={evt => {
                setCurrentFile(!!(evt.currentTarget && evt.currentTarget.files && evt.currentTarget.files[0]))
              }}
              type="file"
            />
            <div
              className={`inline-flex flex-row items-center gap-4 ${!props.output || currentFile ? 'hidden' : ''}`}
              onClick={evt => {
                evt.preventDefault()
                if (fileInputRef.current) {
                  fileInputRef.current.click()
                }
              }}
            >
              <button className="btn btn-lg">Choose File</button>
              <span className="text-lg">{props.output||'No file chosen'}</span>
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
                disabled={!currentFile}
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
  .build()
