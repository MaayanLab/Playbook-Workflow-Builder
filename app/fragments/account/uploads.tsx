import React from 'react'
import { z } from 'zod'
import useSWRMutation from 'swr/mutation'
import { useRouter } from 'next/router'
import dynamic from 'next/dynamic'
import useSWR from 'swr'
import fetcher from '@/utils/next-rest-fetcher'
import * as schema from '@/db'
import type { TypedSchemaRecord } from '@/spec/sql'
import { FileInput, FileURL } from '@/components/core/file'
import { delete_icon, fork_icon } from '@/icons'
import classNames from 'classnames'

const Icon = dynamic(() => import('@/app/components/icon'))
const Bp4Alert = dynamic(() => import('@blueprintjs/core').then(({ Alert }) => Alert))

const deleter = (endpoint: string, { arg }: { arg: any }) => fetch(`${endpoint}/${arg}/delete`, { method: 'POST' }).then(res => res.json())

function humanSize(size: number) {
  const units = ['B', 'KB', 'MB', 'GB', 'TB']
  let unit
  for (let i = 0; i < units.length-1; i++) {
    unit = units[i]
    if (size < 1000) {
      break
    } else {
      size /= 1000
    }
  }
  return `${size.toFixed(2)} ${unit}`
}


export default function Uploads() {
  const router = useRouter()
  const { data: uploads, isLoading } = useSWR<Array<TypedSchemaRecord<typeof schema.user_upload_complete>>>('/api/db/user/uploads', fetcher)
  const [uploadToDelete, setUploadToDelete] = React.useState<TypedSchemaRecord<typeof schema.user_upload_complete> | undefined>(undefined)
  const { trigger: deleteUpload, isMutating } = useSWRMutation('/api/db/user/uploads', deleter)
  return (
    <>
      <h3 className="bp4-heading">Uploads</h3>
      <div className="hero">
        <div className="hero-content text-center">
          <div className="max-w-md">
            <h1 className="text-5xl font-bold prose">Work In Progress</h1>
            <p className="py-6 prose">This feature is still in development. Uploads may not be persistent.</p>
          </div>
        </div>
      </div>
      <progress className={classNames('progress w-full', { 'hidden': !(isLoading || isMutating) })}></progress>
      {uploads ? (
        <div className="overflow-x-auto">
          <table className="table table-compact w-full text-black dark:text-white">
            <thead>
              <tr>
                <th>Filename</th>
                <th>URL</th>
                <th>sha256</th>
                <th>Size</th>
                <th>Timestamp</th>
                <th>Actions</th>
              </tr>
            </thead>
            <tbody>
              {uploads.length === 0 ? <tr><td colSpan={6} align="center">No uploads</td></tr> : null}
              {uploads.map(upload => (
                <tr key={upload.id}>
                  <td>{upload.filename}</td>
                  <td>{upload.url}</td>
                  <td>{upload.sha256.slice(0, 5)}...{upload.sha256.slice(-5)}</td>
                  <td>{humanSize(upload.size)}</td>
                  <td>{upload.created.toString()}</td>
                  <td className="flex flex-row">
                    <button onClick={async () => {
                      const req = await fetch(`/api/db/fpl/start/extend`, {
                        method: 'POST',
                        body: JSON.stringify({
                          type: FileInput.spec,
                          inputs: {},
                          data: {
                            type: FileURL.spec,
                            value: FileURL.codec.encode({ url: upload.url, size: upload.size, filename: upload.filename }),
                          },
                        })
                      })
                      const res = z.string().parse(await req.json())
                      router.push(`/graph/${res}/extend`)
                    }}>
                      <Icon icon={fork_icon} className="fill-black dark:fill-white" />
                    </button>
                    <button onClick={() => {
                      setUploadToDelete(upload)
                    }}>
                      <Icon icon={delete_icon} className="fill-black dark:fill-white" />
                    </button>
                  </td>
                </tr>
              ))}
              <tr><td colSpan={6} align="center">
                <button className="btn btn-primary btn-sm" onClick={async () => {
                  const req = await fetch(`/api/db/fpl/start/extend`, {
                    method: 'POST',
                    body: JSON.stringify({
                      type: FileInput.spec,
                      inputs: {},
                      data: { type: FileURL.spec, },
                    })
                  })
                  const res = z.string().parse(await req.json())
                  router.push(`/graph/${res}`)
                }}>Upload a new file</button>
              </td></tr>
            </tbody>
          </table>
        </div>
      ) : null}
      <Bp4Alert
        cancelButtonText="Cancel"
        confirmButtonText="Delete Upload"
        icon="delete"
        intent="danger"
        isOpen={uploadToDelete !== undefined}
        canEscapeKeyCancel
        canOutsideClickCancel
        onCancel={() => {setUploadToDelete(undefined)}}
        onConfirm={() => {
          if (!uploadToDelete) return
          deleteUpload(uploadToDelete.id, { revalidate: true })
            .then(() => setUploadToDelete(undefined))
        }}
      >
        Are you sure you want to delete {uploadToDelete?.filename} uploaded at {uploadToDelete?.created.toString()}?
        After clicking Delete Upload, your upload will be subject to deletion and <b>cannot be restored</b>.<br />
      </Bp4Alert>
    </>
  )
}