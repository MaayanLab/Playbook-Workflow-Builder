import React from 'react'
import { MetaNode } from '@/spec/metanode'
import { z } from 'zod'

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
      <div>
        <form onSubmit={async (evt) => {
          evt.preventDefault()
          const formData = new FormData(evt.currentTarget)
          const res = await fetch(`/api/components/file/upload`, { method: 'POST', body: formData })
          const records: { file: string[] } = await res.json()
          props.submit(records.file[0])
        }}>
          <input type="file" name="file" />
          <input type="submit" />
        </form>
      </div>
    )
  })
  .build()
