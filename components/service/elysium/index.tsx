import React from 'react'
import { GeneCountMatrix } from '@/components/data/gene_count_matrix'
import { MetaNode } from '@/spec/metanode'
import { z } from 'zod'
import { v4 as uuid4 } from 'uuid'
import python from '@/utils/python'
import SafeRender from '@/utils/saferender'
import * as Auth from 'next-auth/react'
import { useSessionWithId } from '@/app/extensions/next-auth/hooks'
import classNames from 'classnames'

export const FASTQAlignment = MetaNode(`FASTQAlignment`)
  .meta({
    label: 'Align FASTQ Files',
    description: 'Align FASTQ files in the cloud',
  })
  .codec(z.object({
    uid: z.string(),
    filenames: z.string().array(),
    organism: z.string(),
    paired: z.boolean().optional().default(false),
  }))
  .inputs()
  .output(GeneCountMatrix)
  .prompt(props => {
    const { data: session } = useSessionWithId()
    const [uid, _] = React.useState(() => uuid4())
    const organism = 'human'
    const [error, setError] = React.useState<string>()
    return <div>
      {!session || !session.user ? (
        <div className="alert alert-warning shadow-lg block prose max-w-none">
          You are required to &nbsp; <button className="btn btn-sm" onClick={() => {Auth.signIn()}}>sign in</button> &nbsp; to upload files.
        </div>
      ) : (
        <>
          <div className={classNames("alert alert-error shadow-lg block prose max-w-none", { 'hidden': !error })}>
            {error}
          </div>
          <form
            className="inline-flex flex-col"
            onSubmit={(evt) => {
              evt.preventDefault()
              const formData = new FormData(evt.currentTarget)
              const files = formData.getAll('file') as File[]
              const paired = formData.get('paired') === 'on'
              if (paired) {
                const pairs: Record<string, Set<string>> = {}
                files.forEach(file => {
                  const m = /^(.+)_(1|2)\.fastq.gz$/.exec(file.name)
                  if (!m) return
                  const prefix = m[1] as string
                  const suffix = m[2] as string
                  if (!(prefix in pairs)) pairs[prefix] = new Set()
                  pairs[prefix].add(suffix)
                })
                if (Object.keys(pairs).some(prefix => pairs[prefix].size < 2)) {
                  setError('All reads should be paired with two files with suffixes _1 & _2')
                  return
                }
              }
              Promise.all(files.map(async (file) => {
                const charonSearchParams = new URLSearchParams()
                charonSearchParams.append('uid', uid)
                const charonReq = await fetch(`/api/v1/components/service/elysium/charon?${charonSearchParams.toString()}`)
                if (!charonReq.ok) throw new Error(`${await charonReq.text()}`)
                const charonRes = await charonReq.json()
                const charonFormData = new FormData()
                for (const k in charonRes.formData) {
                  charonFormData.append(k, charonRes.formData[k])
                }
                charonFormData.append('file', file)
                // TODO: progress bar (?) https://github.com/MaayanLab/biojupies/blob/6d087aebbd2ff5bc0e8eace5cb67c1f5e19d85e3/website/app/templates/upload/upload_reads.html#L139
                const s3Req = await fetch(charonRes.url, { method: 'POST', body: charonFormData })
                if (!s3Req.ok) throw new Error(`Failed to upload ${file.name}: ${await s3Req.text()}`)
                return file.name
              })).then(filenames => {props.submit({ uid, filenames, paired, organism }, false)})
              .catch(err => {
                console.error(err.toString())
                setError(err.toString())
              })
            }}
          >
            <label className="bp5-control bp5-switch">
              <input type="checkbox" name="paired" />
              <span className="bp5-control-indicator"></span>
              Paired-end reads (pairs should be marked with suffixes: _1 & _2)
            </label>
            <input type="file" name="file" multiple accept=".fastq.gz" />
            <input type="submit" />
          </form>
          {props.output ? <SafeRender component={GeneCountMatrix.view} props={props.output} /> : null}
        </>
      )}
    </div>
  })
  // NOTE: this could also be done in prompt with an API call potentially..
  .resolve(async (props) => await python(
    'components.service.elysium.process_alignment',
    { kargs: [props.data.uid, props.data.filenames], kwargs: { paired: props.data.paired, organism: props.data.organism } },
    message => props.notify({ type: 'info', message }),
  ))
  .story(props => ({
    abstract: `Uploaded FASTQ files were uniformly aligned using Elysium\\ref{doi:10.1101/382937}`,
    legend: `A preview of the aligned gene count matrix.`,
  }))
  .build()
