import dynamic from 'next/dynamic'
import Head from 'next/head'
import Image from 'next/image'
import qrcode from './qrcode.svg'
import landingPageScreenShot from './LandingPageScreenshot.png'

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
          <h5 className="prose"><a className="text-blue-700 dark:text-blue-200 no-underline" href="https://playbook-workflow-builder.cloud">playbook-workflow-builder.cloud</a></h5>
        </div>
        <div className="flex flex-row flex-wrap justify-center mx-auto px-5 my-4 gap-4">
          <div className="flex flex-col items-center gap-2" style={{ flex: '2 0 40em' }}>
            <div className="prose text-xl">
              <h1 className="text-cyan-700 my-0 ml-2">OBJECTIVES</h1>
              <p className="p-0">
                The Playbook Workflow Builder (PWB) is a web-based platform that
                facilitates the interactive construction of workflows by enabling users to
                utilize an ever-growing network of input datasets, semantically annotated
                API endpoints, and data visualization tools contributed by an ecosystem. Via
                a user-friendly web-based user interface (UI), workflows can be constructed
                from contributed building blocks without technical expertise. The workshop
                will demonstrate how to use the platform via interactive sessions, as well as,
                introduce ways the community can contribute to the project.
              </p>
              <Image className="p-2 m-0" src={landingPageScreenShot} alt="Screenshot of PWB Landing Page" objectFit='cover' />
            </div>
            <div className="prose text-xl">
              <h1 className="text-cyan-700 my-0 ml-2">AGENDA</h1>
              <h2 className="my-0">How to use the system to build workflows?</h2>
              <ul>
                <li>Introduction to the Playbook Workflow Builder</li>
                <li>Using the PWB to Create Bioinformatics Workflows</li>
                <li>Published Use Case Workflows</li>
              </ul>
              <h2 className="my-0">How to contribute components and workflows to expand the system?</h2>
              <ul>
                <li>Contributing to the Knowledge Resolution Network of Microservices</li>
                <li>Discussion, demos, interactive sessions, and Q&A</li>
              </ul>
            </div>
          </div>
          <div className="flex flex-col items-center gap-2" style={{ flex: '1 0 20em' }}>
            <div className="card bg-secondary rounded-xl p-4 prose w-96">
              <div className="card-title justify-center text-3xl">
                WHERE
              </div>
              <div className="card-body items-center text-xl">
                Virtual, Requires Registration
              </div>
              <div className="card-title justify-center text-3xl">
                WHEN
              </div>
              <div className="card-body items-center text-xl">
                Tuesday, June 25<br />
                1:00 PM - 4:00 PM ET
              </div>
            </div>
            <div className="card bg-secondary rounded-xl p-4 prose w-96">
              <div className="card-title justify-center text-3xl">
                REGISTRATION
              </div>
              <div className="card-body items-center">
                <Image className="m-0" src={qrcode} alt="Registration Link QR Code" objectFit='cover' />
                <a href="https://mssm.zoom.us/webinar/register/WN_nB40DQScS9y2klZ6frk8JQ#/registration">mssm.zoom.us/webinar/register/WN_<br />nB40DQScS9y2klZ6frk8JQ#/registration</a>
              </div>
            </div>
            <div className="card bg-secondary rounded-xl p-4 prose w-96">
              <div className="card-title justify-center text-3xl">
                FIND OUT MORE
              </div>
              <div className="card-body break-words p-2">
                <h3 className="underline my-0">User Guide</h3>
                <ul className="m-0 py-0"><li className="m-0 p-0"><a href="https://github.com/MaayanLab/Playbook-Workflow-Builder/blob/main/docs/user/index.md">github.com/MaayanLab/Playbook-Workflow-Builder/blob/main/docs/user/index.md</a></li></ul>
              </div>
              <div className="card-body break-words p-2">
                <h3 className="underline my-0">Developer Guide</h3>
                <ul className="m-0 py-0"><li className="m-0 p-0"><a href="https://github.com/MaayanLab/Playbook-Workflow-Builder/blob/main/docs/index.md">github.com/MaayanLab/Playbook-Workflow-Builder/blob/main/docs/index.md</a></li></ul>
              </div>
              <div className="card-body break-words p-2">
                <h3 className="underline my-0">Publication</h3>
                <ul className="m-0 py-0"><li className="m-0 p-0"><a href="https://www.biorxiv.org/content/10.1101/2024.06.08.598037v1">www.biorxiv.org/content/10.1101/2024.06.08.598037v1</a></li></ul>
              </div>
              <div className="card-body break-words p-2">
                <h3 className="underline my-0">Use Cases</h3>
                <ul className="m-0 py-0"><li className="m-0 p-0"><a href="https://www.youtube.com/playlist?list=PLfq4yYrYksVhoPn7xrPQrVYApU-rVLBv1">www.youtube.com/playlist?list=PLfq4yYrYksVhoPn7xrPQrVYApU-rVLBv1</a></li></ul>
              </div>
            </div>
          </div>
        </div>
      </main>
    </>
  )
}
