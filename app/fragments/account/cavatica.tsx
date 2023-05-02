import Image from 'next/image'
import dynamic from 'next/dynamic'
import CAVATICAAPIKeyGuide from '@/app/public/CAVATICA-guide-apikey.png'
import CAVATICAProjectGuide from '@/app/public/CAVATICA-guide-project.png'

const Bp4FormGroup = dynamic(() => import('@blueprintjs/core').then(({ FormGroup }) => FormGroup))
const Bp4InputGroup = dynamic(() => import('@blueprintjs/core').then(({ InputGroup }) => InputGroup))
const Bp4ControlGroup = dynamic(() => import('@blueprintjs/core').then(({ ControlGroup }) => ControlGroup))
const Bp4Button = dynamic(() => import('@blueprintjs/core').then(({ Button }) => Button))

export default function CAVATICA() {
  return (
    <>
      <h3 className="bp4-heading">CAVATICA Integration</h3>
      <div className="hero">
        <div className="hero-content text-center">
          <div className="max-w-md">
            <h1 className="text-5xl font-bold">Coming Soon</h1>
            <p className="py-6 prose">This feature is currently in development. This is currently a non-functioning mockup.</p>
          </div>
        </div>
      </div>
      <p className="bp4-running-text">
        <a href="https://www.cavatica.org/">CAVATICA</a> is a data analysis and sharing platform designed to accelerate discovery in a scalable, cloud-based compute environment where data, results, and workflows are shared among the world's research community. Developed by Seven Bridges and funded in-part by a grant from the National Institutes of Health (NIH) Common Fund, CAVATICA is continuously updated with new tools and datasets.</p>
        <p className="my-2">CAVATICA offers secure storage and compute in a cloud environment. Appyter integration with CAVATICA enables you to execute Appyters against data in a CAVATICA project and use CAVATICA-managed computational resources.</p>
        <p className="my-2">To use CAVATICA, you must <a href="https://pgc-accounts.sbgenomics.com/auth/login" target="_blank">login or create an account</a>, then a project should be established which will be used for file storage and executions that are configured.</p>
        <Image src={CAVATICAProjectGuide} alt="CAVATICA Project Creation Guide" />
        <p className="my-2">To use CAVATICA you must register your CAVATICA API Key with Appyters, this key can be located at <a href="https://cavatica.sbgenomics.com/developer/token" target="_blank">https://cavatica.sbgenomics.com/developer/token</a></p>
        <Image src={CAVATICAAPIKeyGuide} alt="CAVATICA API Key Guide" />
        <div className="my-4">
          <Bp4FormGroup
            label="API Key"
            labelFor="cavatica-api-key"
            helperText="Provide your CAVATICA credentials"
          >
            <Bp4InputGroup
              id="cavatica-api-key"
              type="string"
              placeholder="e.g. 08cd35123..."
              leftIcon="key"
            />
          </Bp4FormGroup>
          <Bp4FormGroup
            label="CAVATICA Project"
            labelFor="cavatica-project"
            helperText="Specify the default project in CAVATICA"
          >
            <Bp4InputGroup
              id="cavatica-project"
              type="string"
              placeholder="e.g. youruser/yourproject"
              leftIcon="projects"
            />
          </Bp4FormGroup>
          <Bp4ControlGroup>
            <Bp4Button
              intent="success"
              onClick={() => {
                // TODO
              }}>Save</Bp4Button>
            <Bp4Button
              intent="danger"
              onClick={() => {
                // TODO
              }}>Delete</Bp4Button>
          </Bp4ControlGroup>
        </div>
    </>
  )
}
