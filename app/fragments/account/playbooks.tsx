import dynamic from 'next/dynamic'
import { delete_icon, edit_icon, fork_icon, save_icon, share_icon, view_report_icon } from '@/icons'

const Icon = dynamic(() => import('@/app/components/icon'))
const Bp4FormGroup = dynamic(() => import('@blueprintjs/core').then(({ FormGroup }) => FormGroup))
const Bp4InputGroup = dynamic(() => import('@blueprintjs/core').then(({ InputGroup }) => InputGroup))

export default function Playbooks() {
  return (
    <>
      <h3 className="bp4-heading">Playbooks</h3>
      <div className="hero">
        <div className="hero-content text-center">
          <div className="max-w-md">
            <h1 className="text-5xl font-bold">Coming Soon</h1>
            <p className="py-6 prose">This feature is currently in development. This is currently a non-functioning mockup.</p>
          </div>
        </div>
      </div>
      <div className="overflow-x-auto">
        <table className="table table-compact w-full">
          <thead>
            <tr>
              <th>Playbook</th>
              <th>Title</th>
              <th>Inputs</th>
              <th>Output</th>
              <th>Timestamp</th>
              <th>Actions</th>
            </tr>
          </thead> 
          <tbody>
            {[
              { id: 'cb557965...', title: 'Use Case 4 Playbook', inputs: 'Gene', output: 'Ranked Tissues', created: '2023-02-12', }, 
              { id: 'b6ed5f65...', title: 'Untitled Playbook', inputs: 'Drug', output: 'Ranked Gene', created: '2023-02-13', }, 
              { id: 'b29c75ec...', title: 'Drugs targeting disease signatures', inputs: 'Disease', output: 'Ranked Drugs', created: '2023-02-14', },
            ].map(playbook => (
              <tr key={playbook.id}>
                <td>{playbook.id}</td>
                <td>{playbook.title}</td>
                <td>{playbook.inputs}</td>
                <td>{playbook.output}</td>
                <td>{playbook.created.toString()}</td>
                <td className="flex flex-row">
                  <button>
                    <Icon icon={share_icon} color="black" />
                  </button>
                  <button>
                    <Icon icon={edit_icon} color="black" />
                  </button>
                  <button>
                    <Icon icon={view_report_icon} color="black" />
                  </button>
                  <button>
                    <Icon icon={fork_icon} color="black" />
                  </button>
                  <button>
                    <Icon icon={delete_icon} color="black" />
                  </button>
                </td>
              </tr>
            ))}
            <tr><td colSpan={6}>&nbsp;</td></tr>
            <tr>
              <td>&nbsp;</td>
              <td>
                <Bp4FormGroup
                  label="Playbook Title"
                  labelInfo="(required)"
                  helperText="A succict title for this playbook"
                >
                  <Bp4InputGroup
                    placeholder="My playbook"
                    leftIcon="label"
                  />
                </Bp4FormGroup>
              </td>
              <td>
                <Bp4FormGroup
                  label="Component Inputs"
                  labelInfo="(required)"
                  helperText="A the inputs to this component"
                >
                  <Bp4InputGroup
                    placeholder="Some component"
                    leftIcon="many-to-one"
                  />
                </Bp4FormGroup>
              </td>
              <td>
                <Bp4FormGroup
                    label="Component Output"
                    labelInfo="(required)"
                    helperText="A the output to this component"
                  >
                    <Bp4InputGroup
                      placeholder="Some component"
                      leftIcon="many-to-one"
                    />
                </Bp4FormGroup>
              </td>
              <td>
                <Bp4FormGroup
                    label="Component Output"
                    labelInfo="(required)"
                    helperText="A the output to this component"
                  >
                    <Bp4InputGroup
                      placeholder="Some component"
                      value={(new Date()).toString()}
                      readOnly
                      leftIcon="many-to-one"
                    />
                </Bp4FormGroup>
              </td>
              <td className="flex flex-row">
                <button>
                  <Icon icon={save_icon} color="black" />
                </button>
                <button>
                  <Icon icon={delete_icon} color="black" />
                </button>
              </td>
            </tr>
          </tbody>
        </table>
      </div>
    </>
  )
}
