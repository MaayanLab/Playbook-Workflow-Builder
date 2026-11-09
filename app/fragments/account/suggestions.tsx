import type { Suggestion } from '@/db'
import { useExRouter } from '@/app/fragments/ex-router'
import React from 'react'
import { useMutation, useQuery, useQueryClient } from '@tanstack/react-query'
import { delete_icon, fork_icon } from '@/icons'
import { useTRPC } from '@/app/trpc/client'
import dynamic from 'next/dynamic'
import { ExLink } from '@/app/fragments/ex-router'
import classNames from 'classnames'

const Icon = dynamic(() => import('@/app/components/icon'))
const Bp5Alert = dynamic(() => import('@blueprintjs/core').then(({ Alert }) => Alert))

export default function Suggestions() {
  const router = useExRouter()
  const trpc = useTRPC()
  const queryClient = useQueryClient()
  const { data: suggestions, isLoading } = useQuery(trpc.users.suggestions.queryOptions())
  const [suggestionToDelete, setSuggestionToDelete] = React.useState<Suggestion | undefined>(undefined)
  const deleteSuggestion = useMutation(trpc.users.deleteSuggestion.mutationOptions({
    onSuccess: () => queryClient.invalidateQueries({ queryKey: trpc.users.suggestions.queryKey() }),
  }))
  const extend = useMutation(trpc.fpprg.extend.mutationOptions())
  return (
    <>
      <h3 className="bp6-heading">Suggestions</h3>
      <progress className={classNames('progress w-full', { 'hidden': !(isLoading || deleteSuggestion.isPending || extend.isPending) })}></progress>
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
                      const res = await extend.mutateAsync({ fpl_id: 'start', process: {
                          type: suggest.name,
                          inputs: {},
                        } })
                      router.push(`/graph/${res}/extend`, undefined, { shallow: true })
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
                <ExLink href="/graph/start/node/start/suggest" shallow><button className="btn btn-primary btn-sm">Suggest a core data type</button></ExLink>
              </td></tr>
            </tbody>
          </table>
        </div>
      ) : null}
      <Bp5Alert
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
          deleteSuggestion.mutateAsync({ id: suggestionToDelete.id })
            .then(() => {
              setSuggestionToDelete(undefined)
            })
        }}
      >
        Are you sure you want to delete {suggestionToDelete?.name} suggestioned at {suggestionToDelete?.created.toString()}?
        After clicking Delete Suggestion, your suggestion will be subject to deletion and <b>cannot be restored</b>.<br />
      </Bp5Alert>
    </>
  )
}
