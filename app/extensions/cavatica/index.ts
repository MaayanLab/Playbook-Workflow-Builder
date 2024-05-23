import uuid5 from "@/utils/uuid"
import packageJson from '@/package.json'

const docker_tag = 'maayanlab/playbook-partnership'
const version = packageJson.version
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
  "s:codeRepository": "https://github.com/MaayanLab/Playbook-Workflow-Builder",
  "s:keywords": ["pwb"],
  "inputs": [
    {
      "id": "config",
      "inputBinding": {
        "position": 1,
        "separate": true,
        "shellQuote": true
      },
      "type": "string",
      "label": "A JSON serialized configuration object"
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

export async function *run_wes_worker({
  url,
  session_id,
  auth_token,
  project,
  api_endpoint = 'https://cavatica-api.sbgenomics.com',
  wes_endpoint = 'wes://cavatica-ga4gh-api.sbgenomics.com',
  polling_interval = 5,
}: {
  url: string,
  session_id: string,
  auth_token: string,
  project: string,
  api_endpoint?: string,
  wes_endpoint?: string,
  polling_interval?: number,
}) {
  const wes_url = wes_endpoint.replace(/^wes:\/\/(.+)$/, 'https://$1/ga4gh/wes')
  const headers = { 'Accept': 'application/json', 'X-SBG-Auth-Token': auth_token }
  // Step 1: Get the right CWL id revision
  const revised_cwl = {...cwl}
  const req0 = await fetch(`${api_endpoint}/v2/apps/${project}/${revised_cwl['id']}`, { headers })
  if (req0.status == 404) revised_cwl['id'] = `${revised_cwl['id']}/0`
  else {
    const res0 = await req0.json()
    if (uuid5(res0['raw']['hints']) == uuid5(revised_cwl['hints'])) {
      revised_cwl['id'] = `${revised_cwl['id']}/${res0['revision']}`
    } else {
      revised_cwl['id'] = `${revised_cwl['id']}/${res0['revision']+1}`
    }
  }
  // Step 2: Submit job to WES
  const body = new FormData()
  body.append('workflow_params', new Blob([JSON.stringify({
    inputs: {
      config: JSON.stringify({
        url,
        session_id,
        auth_token,
        project,
      })
    },
    project
  })], { type: 'application/json' }))
  body.append('workflow_url', new Blob(['#/workflow_attachment/0'], { type: 'text/plain' }))
  body.append('workflow_attachment', new Blob([JSON.stringify(revised_cwl)], { type: 'application/json' }), "pwb-worker.cwl")
  body.append('workflow_type', new Blob(['CWL'], { type: 'text/plain' }))
  body.append('workflow_type_version', new Blob([revised_cwl['cwlVersion']], { type: 'text/plain' }))
  const req1 = await fetch(`${wes_url}/v1/runs`, {
    method: 'POST',
    headers,
    body,
  })
  if (!req1.ok) {
    console.error(`${session_id}: Failed to determine CWL id revision (${req1.status})`)
    console.error(await req1.text())
    throw new Error('Failed to submit WES job')
  }
  const res1 = await req1.json()
  const run_id = res1['run_id']
  yield { run_id, state: null }
  // Step 3: Monitor progress
  let current_state = null
  while (true) {
    await sleep(polling_interval * (0.5 + Math.random()))
    let state: string
    const req2 = await fetch(`${wes_url}/v1/runs/${run_id}/status`, { headers, method: 'GET' })
    if (!req2.ok) {
      if (current_state === 'UNKNOWN') {
        state = 'COMPLETE'
      } else {
        state = 'UNKNOWN'
      }
    } else {
      const res2 = await req2.json()
      state = res2['state']
    }
    if (current_state != state) {
      current_state = state
      yield { run_id, state }
    }
    if (['COMPLETE', 'CANCELED', 'EXECUTOR_ERROR'].includes(current_state)) {
      break
    }
  }
}

export async function abort_wes_worker({
  run_id,
  auth_token,
  wes_endpoint = 'wes://cavatica-ga4gh-api.sbgenomics.com',
}: { run_id: string, auth_token: string, wes_endpoint?: string }) {
  const wes_url = wes_endpoint.replace(/^wes:\/\/(.+)$/, 'https://$1/ga4gh/wes')
  const headers = { 'Accept': 'application/json', 'X-SBG-Auth-Token': auth_token }
  const req = await fetch(`${wes_url}/v1/runs/${run_id}/cancel`, { headers, method: 'POST' })
  return await req.json()
}
