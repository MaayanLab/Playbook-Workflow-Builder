/**
 * API keys let users use the API as themselves (e.g. to upload files, see docs/user/api.md).
 *  A key is only shown when it's created: only a hash of it is stored.
 */
import React from 'react'
import { useMutation, useQuery, useQueryClient } from '@tanstack/react-query'
import classNames from 'classnames'
import { authClient } from '@/app/auth-client'
import usePublicUrl from '@/utils/next-public-url'

const EXPIRATIONS = [
  { label: '30 days', seconds: 30 * 24 * 60 * 60 },
  { label: '90 days', seconds: 90 * 24 * 60 * 60 },
  { label: '1 year', seconds: 365 * 24 * 60 * 60 },
  { label: 'Never', seconds: null },
]

const apiKeysQueryKey = ['account', 'api-keys']

export default function ApiKeys() {
  const queryClient = useQueryClient()
  const publicUrl = usePublicUrl({ absolute: true })
  const { data: apiKeys, isLoading } = useQuery({
    queryKey: apiKeysQueryKey,
    queryFn: async () => {
      const { data, error } = await authClient.apiKey.list()
      if (error) throw new Error(error.message)
      return data.apiKeys
    },
  })
  const [name, setName] = React.useState('')
  const [expiresIn, setExpiresIn] = React.useState<number | null>(EXPIRATIONS[1].seconds)
  // the key just created, shown once
  const [created, setCreated] = React.useState<{ name: string, key: string } | undefined>()
  const [copied, setCopied] = React.useState(false)
  const create = useMutation({
    mutationFn: async () => {
      const { data, error } = await authClient.apiKey.create({ name, expiresIn })
      if (error) throw new Error(error.message)
      return data
    },
    onSuccess: (data) => {
      setCreated({ name: data.name ?? name, key: data.key })
      setCopied(false)
      setName('')
      queryClient.invalidateQueries({ queryKey: apiKeysQueryKey })
    },
  })
  const remove = useMutation({
    mutationFn: async (keyId: string) => {
      const { error } = await authClient.apiKey.delete({ keyId })
      if (error) throw new Error(error.message)
    },
    onSuccess: () => queryClient.invalidateQueries({ queryKey: apiKeysQueryKey }),
  })
  return (
    <div className="overflow-y-auto">
      <h3 className="bp6-heading">API Keys</h3>
      <p className="prose max-w-none">
        API keys let programs use the API as you, for example to upload files, see the <a href={`${publicUrl}/api/openapi.json`}>API description</a> and
        the <a href="https://github.com/MaayanLab/Playbook-Workflow-Builder/blob/main/docs/user/api.md">API guide</a>.
        Send a key in the <code>x-api-key</code> header (or as <code>Authorization: Bearer &lt;key&gt;</code>), and keep it secret: it can do anything you can.
      </p>
      <form
        className="flex flex-row flex-wrap items-end gap-2 my-2"
        onSubmit={evt => { evt.preventDefault(); create.mutate() }}
      >
        <label className="form-control">
          <span className="label-text">Name</span>
          <input className="input input-bordered input-sm" required maxLength={32} placeholder="e.g. my laptop" value={name} onChange={evt => setName(evt.currentTarget.value)} />
        </label>
        <label className="form-control">
          <span className="label-text">Expires</span>
          <select className="select select-bordered select-sm" value={String(expiresIn)} onChange={evt => setExpiresIn(evt.currentTarget.value === 'null' ? null : Number(evt.currentTarget.value))}>
            {EXPIRATIONS.map(({ label, seconds }) => <option key={label} value={String(seconds)}>{label}</option>)}
          </select>
        </label>
        <button type="submit" className="btn btn-primary btn-sm" disabled={!name || create.isPending}>Create API key</button>
      </form>
      {create.error ? <div className="alert alert-error my-2">{create.error.message}</div> : null}
      {created ?
        <div className="alert alert-success my-2 flex flex-col items-start">
          <span>Your new API key <b>{created.name}</b>. Copy it now, it won't be shown again:</span>
          <div className="flex flex-row gap-2 items-center w-full">
            <code className="break-all select-all bg-white text-black p-1 rounded">{created.key}</code>
            <button
              type="button"
              className="btn btn-sm"
              onClick={() => { navigator.clipboard.writeText(created.key).then(() => setCopied(true)) }}
            >{copied ? 'Copied' : 'Copy'}</button>
          </div>
        </div>
        : null}
      <progress className={classNames('progress w-full', { 'hidden': !(isLoading || remove.isPending) })}></progress>
      {apiKeys ?
        <table className="table table-compact w-full text-black dark:text-white">
          <thead>
            <tr><th>Name</th><th>Key</th><th>Created</th><th>Expires</th><th>Last used</th><th></th></tr>
          </thead>
          <tbody>
            {apiKeys.length === 0 ? <tr><td colSpan={6} align="center">No API keys</td></tr> : null}
            {apiKeys.map(apiKey => (
              <tr key={apiKey.id}>
                <td>{apiKey.name}</td>
                <td><code>{apiKey.start}…</code></td>
                <td>{new Date(apiKey.createdAt).toLocaleDateString()}</td>
                <td>{apiKey.expiresAt ? new Date(apiKey.expiresAt).toLocaleDateString() : 'Never'}</td>
                <td>{apiKey.lastRequest ? new Date(apiKey.lastRequest).toLocaleString() : 'Never'}</td>
                <td>
                  <button
                    type="button"
                    className="btn btn-sm btn-error btn-outline"
                    onClick={() => { if (confirm(`Revoke the API key "${apiKey.name}"? Programs using it will stop working.`)) remove.mutate(apiKey.id) }}
                  >Revoke</button>
                </td>
              </tr>
            ))}
          </tbody>
        </table>
        : null}
    </div>
  )
}
