import * as schema from '@/db'
import type { TypedSchemaRecord } from '@/spec/sql'
import { useRouter } from 'next/router'
import React from 'react'
import useSWR from 'swr'
import useSWRMutation from 'swr/mutation'
import { delete_icon, fork_icon } from '@/icons'
import { z } from 'zod'
import fetcher, { fetcherPOST } from '@/utils/next-rest-fetcher'
import dynamic from 'next/dynamic'
import Link from 'next/link'
import classNames from 'classnames'

const Icon = dynamic(() => import('@/app/components/icon'))
const Bp4Alert = dynamic(() => import('@blueprintjs/core').then(({ Alert }) => Alert))

export default function Suggestions() {
  const router = useRouter()
  const { data: suggestions, isLoading } = useSWR<Array<TypedSchemaRecord<typeof schema.suggestion>>>('/api/db/user/suggestions', fetcher)
  const [suggestionToDelete, setSuggestionToDelete] = React.useState<TypedSchemaRecord<typeof schema.suggestion> | undefined>(undefined)
  const { trigger: deleteSuggestion, isMutating } = useSWRMutation(() => suggestionToDelete ? `/api/db/user/suggestions/${suggestionToDelete.id}/delete` : null, fetcherPOST)
  return (
    <>
      <h3 className="bp4-heading">Suggestions</h3>
      <progress className={classNames('progress w-full', { 'hidden': !(isLoading || isMutating) })}></progress>
      {suggestions ? (
        <div className="overflow-x-auto">
          <table className="table table-compact w-full text-black dark:text-white">
            <thead>
              <tr>
                <th></th>
                <th>Title</th>
                <th>Inputs</th>
                <th>Output</th>
                <th>Timestamp</th>
                <th>Actions</th>
                <th></th>
              </tr>
            </thead>
            <tbody>
              {suggestions.length === 0 ? <tr><td colSpan={7} align="center">No suggestions registered</td></tr> : null}
              {suggestions.map(suggest => (
                <tr key={suggest.id}>
                  <td></td>
                  <td>{suggest.name}</td>
                  <td>{suggest.inputs}</td>
                  <td>{suggest.output}</td>
                  <td>{suggest.created.toString()}</td>
                  <td className="flex flex-row">
                    <button onClick={async () => {
                      const req = await fetch(`/api/db/fpl/start/extend`, {
                        headers: {
                          'Content-Type': 'application/json',
                        },
                        method: 'POST',
                        body: JSON.stringify({
                          type: suggest.name,
                          inputs: {},
                        })
                      })
                      const res = z.string().parse(await req.json())
                      router.push(`/graph/${res}/extend`)
                    }}>
                      <Icon icon={fork_icon} className="fill-black dark:fill-white" />
                    </button>
                    <button onClick={() => {
                      setSuggestionToDelete(suggest)
                    }}>
                      <Icon icon={delete_icon} className="fill-black dark:fill-white" />
                    </button>
                  </td>
                  <td></td>
                </tr>
              ))}
              <tr><td colSpan={7} align="center">
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
          deleteSuggestion(undefined, { revalidate: true })
            .then(() => setSuggestionToDelete(undefined))
        }}
      >
        Are you sure you want to delete {suggestionToDelete?.name} suggestioned at {suggestionToDelete?.created.toString()}?
        After clicking Delete Suggestion, your suggestion will be subject to deletion and <b>cannot be restored</b>.<br />
      </Bp4Alert>
    </>
  )
}
