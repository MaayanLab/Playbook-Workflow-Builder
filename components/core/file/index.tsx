import React from 'react'
import { MetaNode } from '@/spec/metanode'
import { z } from 'zod'
import dynamic from 'next/dynamic'
import { file_icon, input_icon } from '@/icons'
import * as Auth from 'next-auth/react'
import { useSessionWithId } from '@/app/extensions/next-auth/hooks'

const BpButton = dynamic(() => import('@blueprintjs/core').then(({ Button }) => Button))

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
    const { data: session } = useSessionWithId()
    return (
      <div>
        {session === null ? (
          <div className="alert alert-warning shadow-lg block">
            You are required to &nbsp; <button className="btn btn-sm" onClick={() => {Auth.signIn()}}>sign in</button> &nbsp; to upload files.
          </div>
        ) : (
          <form
            className="my-2"
            onSubmit={async (evt) => {
              evt.preventDefault()
              const formData = new FormData(evt.currentTarget)
              const res = await fetch(`/api/components/core/file/upload`, { method: 'POST', body: formData })
              const records: { file: string[] } = await res.json()
              props.submit(records.file[0])
            }}
          >
            <div className="form-control flex-row items-center border-solid border-black">
              <input type="file" name="file" className="file-input file-input-lg" />
              <div><BpButton type="submit" large>Submit</BpButton></div>
            </div>
          </form>
        )}
      </div>
    )
  })
  .build()
