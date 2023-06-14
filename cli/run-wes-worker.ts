import uuid5 from "@/utils/uuid"

const docker_tag = 'maayanlab/playbook-partnership'
const version = 'c4'
const cwl = {
  "cwlVersion": "v1.2",
  "class": "CommandLineTool",
  "id": "pwb-worker",
  "baseCommand": ["/app/cli/wes-worker.sh"],
  "hints": [
    {
      "class": "DockerRequirement",
      "dockerPull": `${docker_tag}:${version}`
    }
  ],
  "label": "Playbook Workflow Builder Worker",
  "doc": "Execute playbook workflows using CAVATICA resources.",
  "s:version": version,
  "s:author": "Daniel J. B. Clarke",
  "s:license": "CC-BY-NC-SA-4.0",
  "s:codeRepository": "https://github.com/nih-cfde/playbook-partnership",
  "s:keywords": ["pwb"],
  "inputs": [
    {
      "id": "socket",
      "inputBinding": {
        "position": 1,
        "separate": true,
        "shellQuote": true
      },
      "type": "string",
      "label": "Location to send realtime update stream"
    },
    {
      "id": "session_id",
      "inputBinding": {
        "position": 2,
        "separate": true,
        "shellQuote": true
      },
      "type": "string",
      "label": "Session id"
    }
  ],
  "outputs": [],
  "$namespaces": {
    "s": "https://schema.org"
  },
  "$schemas": [
    "https://schema.org/version/latest/schemaorg-current-http.rdf"
  ]
}

async function sleep(s: number) {
  return await new Promise<void>((resolve, reject) => {setTimeout(() => resolve(), s*1000)})
}

async function main({
  socket,
  session_id,
  auth_token,
  project,
  api_endpoint = 'https://cavatica-api.sbgenomics.com',
  wes_endpoint = 'wes://cavatica-ga4gh-api.sbgenomics.com',
}: { socket: string, session_id: string, auth_token: string, project: string, api_endpoint?: string, wes_endpoint?: string }) {
  const wes_url = wes_endpoint.replace(/^wes:\/\/(.+)$/, 'https://$1/ga4gh/wes')
  const headers = { 'Accept': 'application/json', 'X-SBG-Auth-Token': auth_token }
  // Step 1: Get the right CWL id revision
  const req0 = await fetch(`${api_endpoint}/v2/apps/${project}/${cwl['id']}`, { headers })
  if (req0.status == 404) cwl['id'] = `${cwl['id']}/0`
  else {
    const res0 = await req0.json()
    if (uuid5(res0['raw']['hints']) == uuid5(cwl['hints'])) {
      cwl['id'] = `${cwl['id']}/${res0['revision']}`
    } else {
      cwl['id'] = `${cwl['id']}/${res0['revision']+1}`
    }
  }
  // Step 2: Submit job to WES
  const body = new FormData()
  body.append('workflow_params', new Blob([JSON.stringify({
    inputs: {
      socket,
      session_id,
    },
    project
  })], { type: 'application/json' }))
  body.append('workflow_url', new Blob(['#/workflow_attachment/0'], { type: 'text/plain' }))
  body.append('workflow_attachment', new Blob([JSON.stringify(cwl)], { type: 'application/json' }), "pwb-worker.cwl")
  body.append('workflow_type', new Blob(['CWL'], { type: 'text/plain' }))
  body.append('workflow_type_version', new Blob([cwl['cwlVersion']], { type: 'text/plain' }))
  const req1 = await fetch(`${wes_url}/v1/runs`, {
    method: 'POST',
    headers,
    body,
  })
  const res1 = await req1.json()
  const run_id = res1['run_id']
  console.log(`Started task with run_id=${run_id}`)
  // Step 3: Monitor progress
  let current_state = null
  while (true) {
    await sleep(5 * (0.5 + Math.random()))
    process.stdout.write('.')
    const req2 = await fetch(`${wes_url}/v1/runs/${run_id}/status`, { headers, method: 'GET' })
    const res2 = await req2.json()
    const state = res2['state']
    if (current_state != state) {
      current_state = state
      process.stdout.write('\n')
      console.log(current_state)
    }
    if (['COMPLETE', 'CANCELED', 'EXECUTOR_ERROR'].includes(current_state)) {
      break
    }
  }
}
main({
  auth_token: process.env.CAVATICA_API_KEY as string,
  project: process.env.CAVATICA_PROJECT as string,
  socket: process.argv[2],
  session_id: process.argv[3],
})
