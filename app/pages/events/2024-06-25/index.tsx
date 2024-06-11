import dynamic from 'next/dynamic'
import Head from 'next/head'

const Contributors = dynamic(() => import('@/app/fragments/playbook/contributors'))

export default function WorkshopPage() {
  return (
    <>
      <Head>
        <title>Playbook Workshop: June 25, 2024</title>
      </Head>
      <main className="flex flex-col justify-items-stretch">
        <div className="">
          <Contributors />
        </div>
        <div className="flex flex-col place-items-center bg-primary p-4">
          <div className="flex flex-row items-center">
            <img
              className="p-2 w-24 cursor-pointer dark:stroke-white"
              src={`/PWB-logo.svg`}
            />
            <h1 className="text-black dark:text-white text-4xl font-bold p-2 cursor-pointer">P<span className="text-2xl">laybook</span> W<span className="text-2xl">orkflow</span> B<span className="text-2xl">uilder</span></h1>
          </div>
          <h1 className="prose text-4xl font-bold">VIRTUAL WORKSHOP 2024</h1>
          <h3 className="prose text-2xl">Learn how to build bioinformatics workflows with no technical expertise</h3>
          <h5 className="prose"><a className="text-blue-700 no-underline" href="https://playbook-workflow-builder.cloud">playbook-workflow-builder.cloud</a></h5>
        </div>
        <div className="flex flex-row flex-wrap justify-center mx-auto px-5 my-4 gap-4">
          <div className="flex flex-col">
            <div className="prose">
              <h1 className="text-cyan-700">OBJECTIVES</h1>
              <p>
                The Playbook Workflow Builder (PWB) is a web-based platform that
                facilitates the interactive construction of workflows by enabling users to
                utilize an ever-growing network of input datasets, semantically annotated
                API endpoints, and data visualization tools contributed by an ecosystem. Via
                a user-friendly web-based user interface (UI), workflows can be constructed
                from contributed building blocks without technical expertise. The workshop
                will demonstrate how to use the platform via interactive sessions, as well as,
                introduce ways the community can contribute to the project.
                <img src="" />
              </p>
            </div>
            <div className="prose">
              <h1 className="text-cyan-700">AGENDA</h1>
              <h2>How to use the system to build workflows?</h2>
              <ul>
                <li>Introduction to the Playbook Workflow Builder</li>
                <li>Using the PWB to Create Bioinformatics Workflows</li>
                <li>Published Use Case Workflows</li>
              </ul>
              <h2>How to contribute components and workflows to expand the system?</h2>
              <ul>
                <li>Contributing a Metanode to the Knowledge</li>
                <li>Resolution Network of Microservices</li>
                <li>Discussion, demos, interactive sessions, and Q&A</li>
              </ul>
            </div>
          </div>
          <div className="flex flex-col gap-2">
            <div className="card bg-secondary rounded-xl p-4 prose">
              <div className="card-title justify-center">
                WHERE
              </div>
              <div className="card-body items-center">
                TBA
              </div>
              <div className="card-title justify-center">
                WHEN
              </div>
              <div className="card-body items-center">
                Tuesday, June 25<br />
                1:00p/10:00a - 4:00p/1:00p ET/PT
              </div>
            </div>
            <div className="card bg-secondary rounded-xl p-4 prose">
              <div className="card-title justify-center">
                REGISTRATION
              </div>
              <div className="card-body items-center">
                TBA
                {/* <Image className="m-0" src={qrcode} alt="https://mssm.zoom.us/meeting/register/tZMucuqsrTgrE9eTf_oNlz-trRd_MkEmRALJ" objectFit='cover' />
                <a href="https://mssm.zoom.us/meeting/register/tZMucuqsrTgrE9eTf_oNlz-trRd_MkEmRALJ">Register Here</a> */}
              </div>
            </div>
            <div className="card bg-secondary rounded-xl p-4 prose">
              <div className="card-title justify-center">
                FIND OUT MORE
              </div>
              <div className="card-body">
                <a href="https://github.com/MaayanLab/Playbook-Workflow-Builder/blob/main/docs/user/index.md">User Guide</a>
                <a href="https://github.com/MaayanLab/Playbook-Workflow-Builder/blob/main/docs/index.md">Developer Guide</a>
                <a href="https://www.biorxiv.org/content/10.1101/2024.06.08.598037v1">Publication</a>
                <a href="https://www.youtube.com/playlist?list=PLfq4yYrYksVhoPn7xrPQrVYApU-rVLBv1">Use Cases</a>
              </div>
            </div>
          </div>
        </div>
      </main>
    </>
  )
}
