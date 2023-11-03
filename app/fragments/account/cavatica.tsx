import React from 'react'
import Image from 'next/image'
import dynamic from 'next/dynamic'
import CAVATICAAPIKeyGuide from '@/app/public/CAVATICA-guide-apikey.png'
import CAVATICAProjectGuide from '@/app/public/CAVATICA-guide-project.png'
import classNames from 'classnames'
import { useAPIMutation, useAPIQuery } from '@/core/api/client'
import { UserIntegrationsCAVATICA, UserIntegrationsCAVATICAUpdate } from '@/app/api/client'

const Bp5FormGroup = dynamic(() => import('@blueprintjs/core').then(({ FormGroup }) => FormGroup))
const Bp5InputGroup = dynamic(() => import('@blueprintjs/core').then(({ InputGroup }) => InputGroup))
const Bp5ControlGroup = dynamic(() => import('@blueprintjs/core').then(({ ControlGroup }) => ControlGroup))
const Bp5Button = dynamic(() => import('@blueprintjs/core').then(({ Button }) => Button))

export default function CAVATICA() {
  const { data: userIntegrations, isLoading } = useAPIQuery(UserIntegrationsCAVATICA, {})
  const { trigger: setUserIntegrations, isMutating } = useAPIMutation(UserIntegrationsCAVATICAUpdate, {})
  const [userIntegrationsDraft, setUserIntegrationsDraft] = React.useState({
    cavatica_api_key: '',
    cavatica_default_project: '',
  })
  React.useEffect(() => {
    if (userIntegrations) {
      setUserIntegrationsDraft({
        cavatica_api_key: userIntegrations.cavatica_api_key,
        cavatica_default_project: userIntegrations.cavatica_default_project,
      })
    }
  }, [userIntegrations])
  return (
    <>
      <h3 className="bp5-heading">CAVATICA Integration</h3>
      <form
        className="prose prose-lg max-w-none"
        onSubmit={async (evt) => {
          evt.preventDefault()
          setUserIntegrations({ body: userIntegrationsDraft })
        }}
        method="POST"
      >
        <p><a href="https://www.cavatica.org/">CAVATICA</a> is a data analysis and sharing platform designed to accelerate discovery in a scalable, cloud-based compute environment where data, results, and workflows are shared among the world's research community. Developed by Seven Bridges and funded in-part by a grant from the National Institutes of Health (NIH) Common Fund, CAVATICA is continuously updated with new tools and datasets.</p>
        <p className="my-2">CAVATICA offers secure storage and compute in a cloud environment. Playbook integration with CAVATICA enables you to execute Playbook Worfklows against data in a CAVATICA project and use CAVATICA-managed computational resources.</p>
        <p className="my-2">To use CAVATICA, you must <a href="https://pgc-accounts.sbgenomics.com/auth/login" target="_blank">login or create an account</a>, then a project should be established which will be used for file storage and executions that are configured.</p>
        <Image src={CAVATICAProjectGuide} alt="CAVATICA Project Creation Guide" />
        <Bp5FormGroup
          label="CAVATICA Project"
          labelFor="cavatica-default-project"
          helperText="Specify the default project in CAVATICA"
        >
          <Bp5InputGroup
            id="cavatica-default-project"
            name='cavatica-default-project'
            type="string"
            placeholder="e.g. youruser/yourproject"
            value={userIntegrationsDraft.cavatica_default_project}
            onChange={evt => { setUserIntegrationsDraft(({ cavatica_default_project: _, ...draft }) => ({...draft, cavatica_default_project: evt.target.value })) }}
            leftIcon="projects"
          />
        </Bp5FormGroup>
        <p className="my-2">To use CAVATICA you must register your CAVATICA API Key with The Playbook Workflow Builder, this key can be located at <a href="https://cavatica.sbgenomics.com/developer/token" target="_blank">https://cavatica.sbgenomics.com/developer/token</a></p>
        <Image src={CAVATICAAPIKeyGuide} alt="CAVATICA API Key Guide" />

        <Bp5FormGroup
          label="API Key"
          labelFor="cavatica-api-key"
          helperText="Provide your CAVATICA credentials"
        >
          <Bp5InputGroup
            id="cavatica-api-key"
            type="string"
            name='cavatica-api-key'
            placeholder="e.g. 08cd35123..."
            value={userIntegrationsDraft.cavatica_api_key}
            onChange={evt => { setUserIntegrationsDraft(({ cavatica_api_key: _, ...draft }) => ({...draft, cavatica_api_key: evt.target.value })) }}
            leftIcon="key"
          />
        </Bp5FormGroup>
        <Bp5ControlGroup>
          <Bp5Button
            type="submit"
            intent="success"
          >Save</Bp5Button>
          <Bp5Button
            intent="danger"
            onClick={() => {
              setUserIntegrations({
                body: {
                  cavatica_api_key: '',
                  cavatica_default_project: '',
                }
              }).then(() => {
                setUserIntegrationsDraft({
                  cavatica_api_key: '',
                  cavatica_default_project: '',
                })
              })
            }}>Delete</Bp5Button>
        </Bp5ControlGroup>
        <progress className={classNames('progress w-full', { 'hidden': !(isLoading || isMutating) })}></progress>
      </form>
    </>
  )
}
