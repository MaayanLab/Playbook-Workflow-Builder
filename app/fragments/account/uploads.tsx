import React from 'react'
import { useMutation, useQuery, useQueryClient } from '@tanstack/react-query'
import { useExRouter } from '@/app/fragments/ex-router'
import dynamic from 'next/dynamic'
import type { UserUploadComplete } from '@/db'
import { useTRPC } from '@/app/trpc/client'
import { FileInput, FileURL } from '@/components/core/file'
import { delete_icon, fork_icon } from '@/icons'
import classNames from 'classnames'

const Icon = dynamic(() => import('@/app/components/icon'))
const Bp5Alert = dynamic(() => import('@blueprintjs/core').then(({ Alert }) => Alert))

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
  const router = useExRouter()
  const trpc = useTRPC()
  const queryClient = useQueryClient()
  const { data: uploads, isLoading } = useQuery(trpc.users.uploads.queryOptions())
  const [uploadToDelete, setUploadToDelete] = React.useState<UserUploadComplete | undefined>(undefined)
  const deleteUpload = useMutation(trpc.users.deleteUpload.mutationOptions({
    onSuccess: () => queryClient.invalidateQueries({ queryKey: trpc.users.uploads.queryKey() }),
  }))
  const extend = useMutation(trpc.fpprg.extend.mutationOptions())
  return (
    <>
      <h3 className="bp5-heading">Uploads</h3>
      <progress className={classNames('progress w-full', { 'hidden': !(isLoading || deleteUpload.isPending || extend.isPending) })}></progress>
      {uploads ? (
        <div className="overflow-x-auto">
          <table className="table table-compact w-full text-black dark:text-white">
            <thead>
              <tr>
                <th></th>
                <th>Filename</th>
                <th>URL</th>
                <th>sha256</th>
                <th>Size</th>
                <th>Timestamp</th>
                <th>Actions</th>
                <th></th>
              </tr>
            </thead>
            <tbody>
              {uploads.length === 0 ? <tr><td colSpan={8} align="center">No uploads</td></tr> : null}
              {uploads.map(upload => (
                <tr key={upload.id}>
                  <td></td>
                  <td><a
                    href={`${process.env.NEXT_PUBLIC_URL}/ga4gh/drs/v1/objects/${upload.id}/access/https/data`}
                    download={upload.filename}
                  >{upload.filename}</a></td>
                  <td><a
                    href={`${process.env.NEXT_PUBLIC_URL}/ga4gh/drs/v1/objects/${upload.id}`}
                  >{upload.url}</a></td>
                  <td>{upload.sha256.slice(0, 5)}...{upload.sha256.slice(-5)}</td>
                  <td>{humanSize(upload.size)}</td>
                  <td>{upload.created.toString()}</td>
                  <td className="flex flex-row">
                    <button onClick={async () => {
                      const res = await extend.mutateAsync({ fpl_id: 'start', process: {
                          type: FileInput.spec,
                          inputs: {},
                          data: {
                            type: FileURL.spec,
                            value: FileURL.codec.encode({ url: upload.url, size: upload.size, filename: upload.filename }),
                          },
                        } })
                      router.push(`/graph/${res}/extend`, undefined, { shallow: true })
                    }}>
                      <Icon icon={fork_icon} className="fill-black dark:fill-white" />
                    </button>
                    <button onClick={() => {
                      setUploadToDelete(upload)
                    }}>
                      <Icon icon={delete_icon} className="fill-black dark:fill-white" />
                    </button>
                  </td>
                  <td></td>
                </tr>
              ))}
              <tr><td colSpan={8} align="center">
                <button className="btn btn-primary btn-sm" onClick={async () => {
                  const res = await extend.mutateAsync({ fpl_id: 'start', process: {
                      type: FileInput.spec,
                      inputs: {},
                      data: {
                        type: FileURL.spec,
                        value: null,
                      },
                    } })
                  router.push(`/graph/${res}`, undefined, { shallow: true })
                }}>Upload a new file</button>
              </td></tr>
            </tbody>
          </table>
        </div>
      ) : null}
      <Bp5Alert
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
          deleteUpload.mutateAsync({ id: uploadToDelete.id })
            .then(() => {
              setUploadToDelete(undefined)
            })
        }}
      >
        Are you sure you want to delete {uploadToDelete?.filename} uploaded at {uploadToDelete?.created.toString()}?
        After clicking Delete Upload, your upload will be subject to deletion and <b>cannot be restored</b>.<br />
      </Bp5Alert>
    </>
  )
}
