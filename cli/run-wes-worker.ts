import { run_wes_worker } from '@/app/extensions/cavatica'

async function main() {
  for await (const status of run_wes_worker({
    url: process.argv[2],
    auth_token: process.env.CAVATICA_API_KEY as string,
    project: process.env.CAVATICA_PROJECT as string,
    session_id: process.argv[3],
  })) {
    if (status.state === null) {
      console.log(`Started task with run_id=${status.run_id}`)
    } else {
      console.log(`Task entered state: ${status.state}`)
    }
  }
}
main()
