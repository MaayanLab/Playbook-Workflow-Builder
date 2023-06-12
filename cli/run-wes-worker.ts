import uuid5 from "@/utils/uuid"

const docker_tag = 'maayanlab/pwb-worker'
const version = '0.0.0'
const cwl = {
  "cwlVersion": "v1.2",
  "class": "CommandLineTool",
  "id": "pwb-worker",
  "baseCommand": ["npm","run","wes-worker"],
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
  "s:codeRepository": "https://github.com/MaayanLab/appyter-catalog",
  "s:keywords": ["pwb"],
  "inputs": [
    {
      "id": "socket",
      "inputBinding": {
        "prefix": "",
        "separate": true,
        "shellQuote": true
      },
      "type": "string",
      "label": "Location to send realtime update stream"
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
  auth_token,
  project,
  api_endpoint = 'https://cavatica-api.sbgenomics.com',
  wes_endpoint = 'wes://cavatica-ga4gh-api.sbgenomics.com',
}: { socket: string, auth_token: string, project: string, api_endpoint: string, wes_endpoint: string }) {
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
  body.append('workflow_params', JSON.stringify({ inputs: { socket }, project }), 'application/json')
  body.append('workflow_url', '#/workflow_attachment/0')
  body.append('workflow_attachment', JSON.stringify(cwl), "pwb-worker.cwl")
  body.append('workflow_type', 'CWL')
  body.append('workflow_type_version', cwl['cwlVersion'])
  const req1 = await fetch(`${wes_url}/v1/runs`, {
    method: 'POST',
    headers,
    body,
  })
  const res1 = await req1.json()
  const run_id = res1['run_id']
  // Step 3: Monitor progress
  let current_state = null
  while (true) {
    await sleep(15 * (0.5 + Math.random()))
    const req2 = await fetch(`${wes_url}/v1/runs/${run_id}/status`, { headers, method: 'GET' })
    const res2: any = req2.json()
    const state = res2['state']
    if (current_state != state) {
      current_state = state
      console.log(current_state)
    }
    if (['COMPLETE', 'CANCELED', 'EXECUTOR_ERROR'].includes(current_state)) {
      break
    }
  }
}
main(JSON.parse(process.argv[2]))
