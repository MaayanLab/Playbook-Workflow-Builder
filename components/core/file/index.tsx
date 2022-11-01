import React from 'react'
import { MetaNode } from '@/spec/metanode'
import { z } from 'zod'
import dynamic from 'next/dynamic'

const BpFileInput = dynamic(() => import('@blueprintjs/core').then(({ FileInput }) => FileInput))
const BpButton = dynamic(() => import('@blueprintjs/core').then(({ Button }) => Button))

export const FileURL = MetaNode.createData('FileURL')
  .meta({
    label: 'File URL',
    description: 'An arbitrary file url',
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
  })
  .inputs()
  .output(FileURL)
  .prompt(props => {
    return (
      <form
        onSubmit={async (evt) => {
          evt.preventDefault()
          const formData = new FormData(evt.currentTarget)
          const res = await fetch(`/api/components/core/file/upload`, { method: 'POST', body: formData })
          const records: { file: string[] } = await res.json()
          props.submit(records.file[0])
        }}
      >
        <div><BpFileInput inputProps={{ name: "file" }} /></div>
        <div><BpButton type="submit">Submit</BpButton></div>
      </form>
    )
  })
  .build()
