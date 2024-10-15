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
    const [progress, setProgress] = React.useState<Record<string, string>>({})
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
                const s3Req = await new Promise<{ ok: boolean, text: string }>((resolve, reject) => {
                  const s3Xhr = new XMLHttpRequest()
                  try {
                    s3Xhr.open('POST', charonRes.url, true)
                    s3Xhr.addEventListener('loadstart', evt => {setProgress(progress => ({ ...progress, [file.name]: '0' }))})
                    s3Xhr.addEventListener('loadend', evt => {setProgress(progress => ({ ...progress, [file.name]: '100' }))})
                    s3Xhr.upload.addEventListener('progress', evt => {
                      const pct = Math.ceil((evt.loaded / evt.total) * 100);
                      setProgress(progress => ({ ...progress, [file.name]: `${pct}` }))
                    }, false)
                    s3Xhr.addEventListener('readystatechange', (evt) => {
                      if (s3Xhr.readyState === 4) {
                        resolve({ ok: s3Xhr.status < 300, text: s3Xhr.responseText })
                      }
                    })
                    s3Xhr.send(charonFormData)
                  } catch (e: any) {
                    throw new Error(e.toString())
                  }
                })
                if (!s3Req.ok) throw new Error(`Failed to upload ${file.name}: ${s3Req.text}`)
                return file.name
              })).then(filenames => {props.submit({ uid, filenames, paired, organism }, false)})
              .catch(err => {
                console.error(err.toString())
                setError(err.toString())
              })
            }}
          >
            <input className="file-input" type="file" name="file" multiple accept=".fastq.gz" />
            <div className="form-control">
              <label className="label cursor-pointer">
                <input type="checkbox" name="paired" className="toggle" defaultChecked={props.data?.paired} />
                &nbsp;
                Paired-end reads (pairs should be marked with suffixes: _1 & _2)
              </label>
            </div>
            {Object.keys(progress).map(filename => <div key={filename}><progress className="progress progress-primary w-56" value={+progress[filename]} max={100} />&nbsp;{filename}</div>)}
            <input type="submit" className="btn" title="Submit" />
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
