import handler from '@/utils/next-rest'
import { UnsupportedMethodError } from '@/spec/error'
import packageJson from '@/package.json'

export default handler(async (req, res) => {
  if (req.method !== 'GET') throw new UnsupportedMethodError()
  res.json(JSON.stringify({
    "id": "cloud.playbook-workflow-builder",
    "name": "Playbook Workflow Builder",
    "type": {
      "group": "org.ga4gh",
      "artifact": "drs",
      "version": "1.0.0"
    },
    "description": "This service provides DRS capabilities for the playbook workflow builder.",
    "organization": {
      "name": "Playbook Partnership",
      "url": "https://playbook-workflow-builder.cloud"
    },
    "contactUrl": "mailto:avi.maayan@mssm.edu",
    "documentationUrl": "https://github.com/nih-cfde/playbook-partnership",
    "createdAt": "2023-04-20T12:20:15Z",
    "updatedAt": "2023-09-05T12:20:15Z",
    "environment": "prod",
    "version": packageJson.version
  }))
})
