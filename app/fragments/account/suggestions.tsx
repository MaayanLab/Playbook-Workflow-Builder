import * as schema from '@/db'
import type { TypedSchemaRecord } from '@/spec/sql'
import { useRouter } from 'next/router'
import React from 'react'
import useSWR from 'swr'
import useSWRMutation from 'swr/mutation'
import { delete_icon, fork_icon } from '@/icons'
import { z } from 'zod'
import fetcher from '@/utils/next-rest-fetcher'
import dynamic from 'next/dynamic'
import Link from 'next/link'

const Icon = dynamic(() => import('@/app/components/icon'))
const Bp4Alert = dynamic(() => import('@blueprintjs/core').then(({ Alert }) => Alert))

const deleter = (endpoint: string, { arg }: { arg: any }) => fetch(`${endpoint}/${arg}/delete`, { method: 'POST' }).then(res => res.json())

export default function Suggestions() {
  const router = useRouter()
  const { data: suggestions, isLoading } = useSWR<Array<TypedSchemaRecord<typeof schema.suggestion>>>('/api/db/user/suggestions', fetcher)
  const [suggestionToDelete, setSuggestionToDelete] = React.useState<TypedSchemaRecord<typeof schema.suggestion> | undefined>(undefined)
  const { trigger: deleteSuggestion, isMutating } = useSWRMutation('/api/db/user/suggestions', deleter)
  return (
    <>
      <h3 className="bp4-heading">Suggestions</h3>
      <progress className={`progress w-full ${isLoading || isMutating ? '' : 'hidden'}`}></progress>
      {suggestions ? (
        <div className="overflow-x-auto">
          <table className="table table-compact w-full">
            <thead>
              <tr>
                <th>Title</th>
                <th>Inputs</th> 
                <th>Output</th>
                <th>Timestamp</th>
                <th>Actions</th>
              </tr>
            </thead> 
            <tbody>
              {suggestions.length === 0 ? <tr><td colSpan={5} align="center">No suggestions registered</td></tr> : null}
              {suggestions.map(suggest => (
                <tr key={suggest.id}>
                  <td>{suggest.name}</td>
                  <td>{suggest.inputs}</td>
                  <td>{suggest.output}</td>
                  <td>{suggest.created.toString()}</td>
                  <td className="flex flex-row">
                    <button onClick={async () => {
                      const req = await fetch(`/api/db/fpl/start/extend`, {
                        method: 'POST',
                        body: JSON.stringify({
                          type: suggest.name,
                          inputs: {},
                        })
                      })
                      const res = z.string().parse(await req.json())
                      router.push(`/graph/${res}/extend`)
                    }}>
                      <Icon icon={fork_icon} color="black" />
                    </button>
                    <button onClick={() => {
                      setSuggestionToDelete(suggest)
                    }}>
                      <Icon icon={delete_icon} color="black" />
                    </button>
                  </td>
                </tr>
              ))}
              <tr><td colSpan={5} align="center">
                <Link href="/graph/start/node/start/suggest"><button className="btn btn-primary btn-sm">Suggest a core data type</button></Link>
              </td></tr>
            </tbody>
          </table>
        </div>
      ) : null}
      <Bp4Alert
        cancelButtonText="Cancel"
        confirmButtonText="Delete suggestion"
        icon="delete"
        intent="danger"
        isOpen={suggestionToDelete !== undefined}
        canEscapeKeyCancel
        canOutsideClickCancel
        onCancel={() => {setSuggestionToDelete(undefined)}}
        onConfirm={() => {
          if (!suggestionToDelete) return
          deleteSuggestion(suggestionToDelete.id, { revalidate: true })
            .then(() => setSuggestionToDelete(undefined))
        }}
      >
        Are you sure you want to delete {suggestionToDelete?.name} suggestioned at {suggestionToDelete?.created.toString()}?
        After clicking Delete Suggestion, your suggestion will be subject to deletion and <b>cannot be restored</b>.<br />
      </Bp4Alert>
    </>
  )
}
