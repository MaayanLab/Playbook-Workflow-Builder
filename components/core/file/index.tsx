import React from 'react'
import { MetaNode } from '@/spec/metanode'
import { z } from 'zod'
import dynamic from 'next/dynamic'
import { file_icon, input_icon } from '@/icons'
import * as Auth from 'next-auth/react'
import { useSessionWithId } from '@/app/extensions/next-auth/hooks'
import classNames from 'classnames'
import { clientUploadFile } from  '@/components/core/file/api/upload/client'

const Bp4Button = dynamic(() => import('@blueprintjs/core').then(({ Button }) => Button))

export const FileURL = MetaNode('FileURL')
  .meta({
    label: 'File URL',
    description: 'URL to a File',
    icon: [file_icon],
  })
  .codec(z.object({
    description: z.string().optional(),
    url: z.string(),
    filename: z.string(),
    size: z.number(),
  }))
  .view(({ url }) => (
    <div>
      <h2>File: {url}</h2>
    </div>
  ))
  .build()

export const FileInput = MetaNode('FileInput')
  .meta({
    label: 'Input File',
    description: 'Upload a Data File',
    icon: [input_icon],
  })
  .inputs()
  .output(FileURL)
  .prompt(props => {
    const [currentFile, setCurrentFile] = React.useState<{ description?: string, url?: string, filename?: string, size?: number }>({})
    const fileInputRef = React.useRef<HTMLInputElement>(null)
    const { data: session } = useSessionWithId()
    const output = React.useMemo(() => props.output || { description: undefined, url: undefined, filename: undefined, size: undefined }, [props.output])
    React.useEffect(() => {
      setCurrentFile(output)
    }, [output])
    return (
      <div>
        {!session || !session.user ? (
          <div className="alert alert-warning shadow-lg block">
            You are required to &nbsp; <button className="btn btn-sm" onClick={() => {Auth.signIn()}}>sign in</button> &nbsp; to upload files.
          </div>
        ) : (
          <form
            className="my-2 inline-flex flex-col"
            onSubmit={async (evt) => {
              evt.preventDefault()
              if (currentFile.url && currentFile.url === props.output?.url && currentFile.filename && currentFile.size) {
                props.submit({
                  description: currentFile.description,
                  url: currentFile.url,
                  filename: currentFile.filename,
                  size: currentFile.size,
                })
              } else {
                const formData = new FormData(evt.currentTarget)
                const rawDescription = formData.get('description')
                const description = rawDescription === null ? undefined : rawDescription.toString()
                const { file: [record] } = await clientUploadFile(formData)
                props.submit({ description, ...record })
              }
            }}
          >
            <input
              ref={fileInputRef}
              name="file"
              className={classNames('file-input file-input-lg', { 'hidden': currentFile.url })}
              onChange={evt => {
                const url = evt.currentTarget && evt.currentTarget.files && evt.currentTarget.files[0] ? evt.currentTarget.files[0].name : output.url
                setCurrentFile(({ description, url: _ }) => ({ description, url }))
              }}
              type="file"
            />
            <div
              className={classNames('inline-flex flex-row items-center gap-4', { 'hidden': !currentFile.url })}
              onClick={evt => {
                evt.preventDefault()
                if (fileInputRef.current) {
                  fileInputRef.current.click()
                }
              }}
            >
              <button className="btn btn-lg">Choose File</button>
              <span className="text-lg">{currentFile.filename||currentFile.url||'No file chosen'}</span>
            </div>
            <div className="bp4-input-group">
              <input
                type="text"
                name="description"
                className="bp4-input"
                placeholder={`File description`}
                onChange={evt => {setCurrentFile(({ description: _, url }) => ({ url, description: evt.target.value }))}}
                value={currentFile.description||''}
              />
            </div>
            <div className="inline-flex flex-row">
              <a
                href="/api/components/core/file/example_matrix.tsv"
                download="example_matrix.tsv"
              >
                <Bp4Button
                  large
                  text="Example"
                  rightIcon="bring-data"
                />
              </a>
              <Bp4Button
                large
                disabled={
                  currentFile.url === undefined
                  || (
                    currentFile.description === output.description
                    && currentFile.url === output.url
                  )
                }
                type="submit"
                text="Submit"
                rightIcon="send-to-graph"
              />
            </div>
          </form>
        )}
      </div>
    )
  })
  .story(props =>
    `A file${props.output?.description ? ` containing ${props.output?.description}` : ''} was first uploaded.`
  )
  .build()


  export const MultiFileInputForCTD = MetaNode('MultiFileInputForCTD')
  .meta({
    label: 'CTD - Precalculations',
    description: 'Upload Multiple Data Files for CTD',
    icon: [input_icon],
  })
  .inputs()
  .output(FileURL)
  .prompt(props => { 
    const { data: session } = useSessionWithId()
    const output = React.useMemo(() => props.output || { description: undefined, url: undefined, filename: undefined, size: undefined }, [props.output])

    const [geneListFile, setGeneListFile] = React.useState<{ description?: string, url?: string, filename?: string, size?: number }>({})
    const geneListInputRef = React.useRef<HTMLInputElement>(null)
    React.useEffect(() => {
      setGeneListFile(output)
    }, [output])

    const [adjMatrixFile, setAdjMatrixFile] = React.useState<{ description?: string, url?: string, filename?: string, size?: number }>({})
    const adjMatrixInputRef = React.useRef<HTMLInputElement>(null)
    React.useEffect(() => {
      setAdjMatrixFile(output)
    }, [output])

    return (
      <div>
        {!session || !session.user ? (
          <div className="alert alert-warning shadow-lg block">
            You are required to &nbsp; <button className="btn btn-sm" onClick={() => {Auth.signIn()}}>sign in</button> &nbsp; to upload files.
          </div>
        ) : (
          <form
            className="my-2 inline-flex flex-col"
            onSubmit={async (evt) => {
              evt.preventDefault()
              if (geneListFile.url && geneListFile.url === props.output?.url && geneListFile.filename && geneListFile.size &&
                        adjMatrixFile.url && adjMatrixFile.url === props.output?.url && adjMatrixFile.filename && adjMatrixFile.size) {                
                props.submit({
                  description: geneListFile.description,
                  url: geneListFile.url,
                  filename: geneListFile.filename,
                  size: geneListFile.size,
                })
              } else {
                const formData = new FormData(evt.currentTarget)

                
                for (const pair of formData.entries()) {
                  console.log(`${pair[0]}, ${pair[1]}`);
                }
                /*
                let rawDescription = formData.get('description');
                let description = rawDescription === null ? undefined : rawDescription.toString();
                let { file: [record] } = await clientUploadFile(formData);
                props.submit({ description, ...record });*/
              }
            }}
          >
            <input
              ref={geneListInputRef}
              name="geneSetFile"
              className={classNames('file-input file-input-lg', { 'hidden': geneListFile.url })}
              onChange={evt => {
                const url = evt.currentTarget && evt.currentTarget.files && evt.currentTarget.files[0] ? evt.currentTarget.files[0].name : output.url
                setGeneListFile(({ description, url: _ }) => ({ description, url }))
              }}
              type="file"
              accept=".csv"
            />
            <div
              className={classNames('inline-flex flex-row items-center gap-4', { 'hidden': !geneListFile.url })}
              onClick={evt => {
                evt.preventDefault()
                if (geneListInputRef.current) {
                  geneListInputRef.current.click()
                }
              }}
            >
              <button className="btn btn-lg">Choose Gene list File</button>
              <span className="text-lg">{geneListFile.filename||geneListFile.url||'No file chosen'}</span>
            </div>
            <div className="bp4-input-group">
              <input
                type="text"
                name="description"
                className="bp4-input"
                placeholder={`File description`}
                onChange={evt => {setGeneListFile(({ description: _, url }) => ({ url, description: evt.target.value }))}}
                value={geneListFile.description||''}
              />
            </div>
            <div className="inline-flex flex-row">
              <a
                href="/api/components/core/file/example_geneSetCTD.csv"
                download="example_geneSetCTD.csv"
              >
                <Bp4Button
                    large
                    text="Example"
                    rightIcon="bring-data"
                />
              </a>
            </div>


      
            <input
              ref={adjMatrixInputRef}
              name="adjMatrixFile"
              className={classNames('file-input file-input-lg', { 'hidden': adjMatrixFile.url })}
              onChange={evt => {
                const url = evt.currentTarget && evt.currentTarget.files && evt.currentTarget.files[0] ? evt.currentTarget.files[0].name : output.url
                setAdjMatrixFile(({ description, url: _ }) => ({ description, url }))
              }}
              type="file"
              accept=".csv"
            />
            <div
              className={classNames('inline-flex flex-row items-center gap-4', { 'hidden': !adjMatrixFile.url })}
              onClick={evt => {
                evt.preventDefault()
                if (adjMatrixInputRef.current) {
                  adjMatrixInputRef.current.click()
                }
              }}
            >
              <button className="btn btn-lg">Choose Adj. Matrix File</button>
              <span className="text-lg">{adjMatrixFile.filename||adjMatrixFile.url||'No file chosen'}</span>
            </div>
            <div className="bp4-input-group">
              <input
                type="text"
                name="description"
                className="bp4-input"
                placeholder={`File description`}
                onChange={evt => {setAdjMatrixFile(({ description: _, url }) => ({ url, description: evt.target.value }))}}
                value={adjMatrixFile.description||''}
              />
            </div>
            <div className="inline-flex flex-row">
              <a
                href="/api/components/core/file/example_geneSetCTD.csv"
                download="example_geneSetCTD.csv"
              >
                <Bp4Button
                    large
                    text="Example"
                    rightIcon="bring-data"
                />
              </a>
            </div>





            <Bp4Button
                large
                disabled={
                  (geneListFile.url === undefined || (geneListFile.description === output.description && geneListFile.url === output.url)) ||
                  (adjMatrixFile.url === undefined || (adjMatrixFile.description === output.description && adjMatrixFile.url === output.url))
                }
                type="submit"
                text="Submit"
                rightIcon="send-to-graph"
            />
          </form>
        )}
      </div> 
    )
  })
  .story(props =>
    `A file${props.output?.description ? ` containing ${props.output?.description}` : ''} was first uploaded.`
  )
  .build()