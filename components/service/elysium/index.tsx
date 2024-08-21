import React from 'react'
import { GeneCountMatrix } from '@/components/data/gene_count_matrix'
import { MetaNode } from '@/spec/metanode'
import { z } from 'zod'
import { v4 as uuid4 } from 'uuid'
import python from '@/utils/python'

export const FASTQAlignment = MetaNode(`FASTQAlignment`)
  .meta({
    label: 'Align FASTQ Files',
    description: 'Align FASTQ files in the cloud',
  })
  .codec(z.object({
    uid: z.string(),
    filenames: z.string().array(),
    organism: z.string(),
  }))
  .inputs()
  .output(GeneCountMatrix)
  .prompt(props => {
    const [uid, _] = React.useState(() => uuid4())
    const organism = 'human'
    const [error, setError] = React.useState<string>()
    return <>
      <form onSubmit={(evt) => {
        evt.preventDefault()
        const formData = new FormData(evt.currentTarget)
        const files = formData.getAll('file') as File[]
        Promise.all(files.map(async (file) => {
          const charonSearchParams = new URLSearchParams()
          charonSearchParams.append('uid', uid)
          const charonReq = await fetch(`/api/v1/components/service/elysium/charon?${charonSearchParams.toString()}`)
          const charonRes = await charonReq.json()
          for (const k in charonRes.formData) {
            formData.append(k, charonRes.formData[k])
          }
          // TODO: progress bar (?) https://github.com/MaayanLab/biojupies/blob/6d087aebbd2ff5bc0e8eace5cb67c1f5e19d85e3/website/app/templates/upload/upload_reads.html#L139
          const s3Req = await fetch(charonRes.url, { method: 'POST', body: formData })
          if (!s3Req.ok) throw new Error(`Failed to upload ${file.name}`)
          return file.name
        })).then(filenames => {props.submit({ uid, filenames, organism }, true)})
        .catch(err => setError(err))
      }}>
        <input type="file" name="file" multiple />
        <input type="submit" />
      </form>
    </>
  })
  // NOTE: this could also be done in prompt with an API call potentially..
  .resolve(async (props) => await python(
    'components.service.elysium.process_single_end',
    { kargs: [props.data.uid, props.data.filenames, props.data.organism] },
    message => props.notify({ type: 'info', message }),
  ))
  .story(props => ({
    abstract: `Uploaded FASTQ files were uniformly aligned using Elysium\\ref{doi:10.1101/382937}`,
    legend: `A preview of the aligned gene count matrix.`,
  }))
  .build()
