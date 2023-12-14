import { UserIntegrationsBioComputeAuth } from "@/app/api/client"
import { useAPIQuery } from "@/core/api/client"
import classNames from "classnames"
import { signIn } from "next-auth/react"
import { useRouter } from "next/router"
import React from "react"

export default function BioCompute() {
  const router = useRouter()
  const { data: biocomputeAuth } = useAPIQuery(UserIntegrationsBioComputeAuth, {})
  React.useEffect(() => {
    if (biocomputeAuth?.orcid && biocomputeAuth.biocompute) {
      if (router.query.callback) {
        router.push(router.query.callback as string)
      }
    }
  }, [biocomputeAuth, router.query.callback])
  return (
    <div className="prose">
      <h3 className="bp5-heading">BioCompute Integration</h3>
      <p>
        BioCompute is a standardized way of representing a computational workflow.
        The <a href="https://biocomputeobject.org/">BioCompute portal</a> lets users create, edit, and publish BioCompute Objects (BCOs).
        To publish playbook partnership workflows on the BioCompute portal, you must have an account with the BioCompute portal.
        Both your playbook workflow builder account AND your biocompute account should be associated with a common <a href="https://orcid.org/">ORCID</a>.
      </p>
      {biocomputeAuth ?
        biocomputeAuth.orcid && biocomputeAuth.biocompute ? (
          <div className="alert alert-success my-2">BioCompute Integration is active</div>
        ) : (
          <div className="flex flex-row gap-2 justify-around">
            <a onClick={() => {signIn('orcid')}} className={classNames('btn', { 'btn-disabled': biocomputeAuth.orcid })}>Link Playbook Workflow Builder Account with ORCID</a>
            <a href="https://biocomputeobject.org/profile" target="_blank" className={classNames('btn', { 'btn-disabled': biocomputeAuth.biocompute })}>Create BCO Account/Link ORCID</a>
          </div>
        ) : null}
    </div>
  )
}