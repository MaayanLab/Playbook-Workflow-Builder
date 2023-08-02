import dynamic from 'next/dynamic'
import { delete_icon, view_in_graph_icon, view_report_icon } from '@/icons'
import { useAPIQuery, useAPIMutation } from '@/core/api/client'
import { DeleteUserPlaybook, UserPlaybooks } from '@/app/api/client'
import Link from 'next/link'

const Icon = dynamic(() => import('@/app/components/icon'))

export default function Playbooks() {
  const { data: playbooks, mutate: mutatePlaybooks } = useAPIQuery(UserPlaybooks, {})
  const { trigger: deleteUserPlaybook } = useAPIMutation(DeleteUserPlaybook, {})
  return (
    <>
      <h3 className="bp4-heading">Playbooks</h3>
      <div className="overflow-x-auto">
        <table className="table table-compact w-full text-black dark:text-white">
          <thead>
            <tr>
              <th></th>
              <th>Title</th>
              <th>Inputs</th>
              <th>Outputs</th>
              <th>Created</th>
              <th>Public</th>
              <th>Clicks</th>
              <th>Actions</th>
              <th></th>
            </tr>
          </thead>
          <tbody>
            {(playbooks ?? []).length === 0 ? <tr><td colSpan={9} align="center">No playbooks saved</td></tr> : null}
            {(playbooks ?? []).map(playbook => (
              <tr key={playbook.id}>
                <td></td>
                <td>{playbook.title}</td>
                <td>{playbook.inputs}</td>
                <td>{playbook.outputs}</td>
                <td>{playbook.created.toString()}</td>
                <td>{playbook.public ? 'Yes' : 'No'}</td>
                <td>{playbook.clicks}</td>
                <td className="flex flex-row gap-2">
                  <Link href={`/report/${playbook.playbook}`}>
                    <button>
                      <Icon icon={view_report_icon} className="fill-black dark:fill-white" />
                    </button>
                  </Link>
                  <Link href={`/graph/${playbook.playbook}`}>
                    <button>
                      <Icon icon={view_in_graph_icon} className="fill-black dark:fill-white" />
                    </button>
                  </Link>
                  <button onClick={() => {
                    deleteUserPlaybook({ query: { id: playbook.playbook }, body: {} })
                      .then(() => {
                        mutatePlaybooks(data => data ? data.filter(({ id }) => id !== playbook.id) : data)
                      })
                  }}>
                    <Icon icon={delete_icon} className="fill-black dark:fill-white" />
                  </button>
                </td>
                <td></td>
              </tr>
            ))}
          </tbody>
        </table>
      </div>
    </>
  )
}
