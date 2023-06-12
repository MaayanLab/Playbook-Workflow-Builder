import re
import io
import json
import time
import click
import random
import requests

docker_tag = 'maayanlab/pwb-worker'
version = '0.0.0'

cwl = {
  "cwlVersion": "v1.2",
  "class": "CommandLineTool",
  "id": "pwb-worker",
  "baseCommand": ["npm","run","wes-worker"],
  "hints": [
    {
      "class": "DockerRequirement",
      "dockerPull": f"{docker_tag}:{version}"
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
        "separate": True,
        "shellQuote": True
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

def ValidateJSON(ctx, param, value):
  if isinstance(value, dict):
    return value
  try:
    return json.loads(value)
  except:
    raise click.BadParameter('Expected json')

@click.command()
@click.option('--inputs', type=click.UNPROCESSED, callback=ValidateJSON, required=True)
@click.option('--auth-token', type=str, required=True)
@click.option('--project', type=str, required=True)
@click.option('--api-endpoint', type=str, default='https://cavatica-api.sbgenomics.com', required=True)
@click.option('--wes-endpoint', type=str, default='wes://cavatica-ga4gh-api.sbgenomics.com', required=True)
def run_in_cavatica(
  inputs: dict,
  auth_token: str,
  project: str,
  api_endpoint: str,
  wes_endpoint: str,
):
  wes_url = re.sub(r'^wes://(.+)$', r'https://$1/ga4gh/wes', wes_endpoint)
  headers=dict({ 'Accept': 'application/json', 'X-SBG-Auth-Token': auth_token })
  # Step 1: Get the right CWL id revision
  req = requests.get(f"{api_endpoint}/v2/apps/{project}/{cwl['id']}", headers=headers)
  if req.status == 404:
    cwl['id'] = f"{cwl['id']}/0"
  else:
    res = req.json()
    if json.dumps(res['raw']['hints'], sort_keys=True) == json.dumps(cwl['hints'], sort_keys=True):
      cwl['id'] = f"{cwl['id']}/{res['revision']}"
    else:
      cwl['id'] = f"{cwl['id']}/{res['revision']+1}"
  # Step 2: Submit job to WES
  req = requests.post(f"{wes_url}/v1/runs", headers=headers, files=dict(
    workflow_params=(None, json.dumps(dict(inputs=inputs, project=project)), 'application/json'),
    workflow_url=(None, '#/workflow_attachment/0'),
    workflow_attachment=("pwb-worker.cwl", io.BytesIO(json.dumps(cwl).encode())),
    workflow_type=(None, 'CWL'),
    workflow_type_version=(None, cwl['cwlVersion']),
  ))
  res = req.json()
  run_id = res['run_id']
  # Step 3: Monitor progress
  current_state = None
  while True:
    time.sleep(15 * (0.5 + random.random()))
    req = requests.get(f"{wes_url}/v1/runs/{run_id}/status", headers=headers)
    res = req.json()
    state = res['state']
    if current_state != state:
      current_state = state
      yield current_state
    if current_state in {'COMPLETE', 'CANCELED', 'EXECUTOR_ERROR'}:
      break

if __name__ == '__main__':
  run_in_cavatica()
