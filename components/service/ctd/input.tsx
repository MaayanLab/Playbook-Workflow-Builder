import React from 'react'
import { MetaNode } from '@/spec/metanode'
import { z } from 'zod'
import dynamic from 'next/dynamic'
import { file_icon, input_icon } from '@/icons'
import * as Auth from 'next-auth/react'
import { useSessionWithId } from '@/app/extensions/next-auth/hooks'
import classNames from 'classnames'
import { clientUploadFile } from  '@/components/core/file/api/upload/client'
import { Table, Cell, Column} from '@/app/components/Table'
import { FileC } from  '@/components/core/file/'

const Bp4Button = dynamic(() => import('@blueprintjs/core').then(({ Button }) => Button))

export const CTDPrecalculationsFileURLs = MetaNode('CTDPrecalculationsFileURLs')
  .meta({
    label: 'CTDPrecalculationsFileURLs',
    description: 'URLs to the files for CTD Precalculations API call',
    icon: [file_icon],
  })
  .codec(z.object({
    geneListFile: FileC,
    adjMatrixFile: FileC,
  }))
  .view(fileURLs => {
    const filenames = ['geneListFile', 'adjMatrixFile'] as const
    return (
      <Table
      height={500}
      cellRendererDependencies={[fileURLs]}
      numRows={filenames.length}
      enableGhostCells
      enableFocusedCell
      >
      <Column
        name="File Name"
        cellRenderer={row => <Cell key={row+''}>{fileURLs[filenames[row]].filename}</Cell>}
      />
      <Column
        name="URL"
        cellRenderer={row => <Cell key={row+''}>{fileURLs[filenames[row]].url}</Cell>}
      />
      <Column
        name="Size"
        cellRenderer={row => <Cell key={row+''}>{fileURLs[filenames[row]].size}</Cell>}
      />
      <Column
        name="Description"
        cellRenderer={row => <Cell key={row+''}>{fileURLs[filenames[row]].description}</Cell>}
      />
      </Table>
    )
  })
  .build()

export const CTDPrecalculationsFileInput = MetaNode('CTDPrecalculationsFileInput')
  .meta({
    label: 'CTD - Precalculations',
    description: 'Upload two files, Gene list and Adj. Matrix, for CTD Precalculations',
    icon: [input_icon],
  })
  .inputs()
  .output(CTDPrecalculationsFileURLs)
  .prompt(props => { 
    const { data: session } = useSessionWithId()

    const [geneListFile, setGeneListFile] = React.useState<z.infer<typeof FileC>>({ url: '', filename: '' })
    const geneListInputRef = React.useRef<HTMLInputElement>(null)

    const [adjMatrixFile, setAdjMatrixFile] = React.useState<z.infer<typeof FileC>>({ url: '', filename: '' })
    const adjMatrixInputRef = React.useRef<HTMLInputElement>(null)

    React.useEffect(() => {
      if (!props.output) return
      setGeneListFile(props.output.geneListFile)
      setAdjMatrixFile(props.output.adjMatrixFile)
    }, [props.output])
    return (
      <div>
        {!session || !session.user ? (
          <div className="alert alert-warning shadow-lg block prose">
            You are required to &nbsp; <button className="btn btn-sm" onClick={() => {Auth.signIn()}}>sign in</button> &nbsp; to upload files.
          </div>
        ) : (
          <form
            className="my-2 inline-flex flex-col"
            onSubmit={async (evt) => {
              evt.preventDefault();
              if (props.output != null && geneListFile.url && geneListFile.url === props.output.geneListFile?.url && geneListFile.filename && geneListFile.size &&
                        adjMatrixFile.url && adjMatrixFile.url === props.output.adjMatrixFile?.url && adjMatrixFile.filename && adjMatrixFile.size) {                
                props.submit({ geneListFile, adjMatrixFile })
              }else{
                const formData = new FormData(evt.currentTarget)
                formData.delete("description");
                /*
                let rawDescription = formData.getAll('description').toString();
                formData.delete("description");
                formData.set("description", rawDescription)
                let description = rawDescription === null ? undefined : rawDescription.toString(); */ 
                let response = await clientUploadFile(formData);   
                const geneSetFileResponse = response.geneSetFile[0];
                const adjMatrixFileResponse = response.adjMatrixFile[0];             

                props.submit({ geneListFile: {
                    description: "Gene list file",
                    url: geneSetFileResponse.url,
                    filename: geneSetFileResponse.filename,
                    size: geneSetFileResponse.size,
                  }, adjMatrixFile: {
                    description: "Adj. Matrix File",
                    url: adjMatrixFileResponse.url,
                    filename: adjMatrixFileResponse.filename,
                    size: adjMatrixFileResponse.size,
                  }
                })
              }
            }}
          >
            <p style={{ fontSize: "16px" }}><b>Select gene list file</b>, acceptable file types: <b>.csv, .txt</b></p>
            <input
              ref={geneListInputRef}
              name="geneSetFile"
              className={classNames('file-input file-input-lg', { 'hidden': geneListFile.url })}
              onChange={evt => {
                const url = evt.currentTarget && evt.currentTarget.files && evt.currentTarget.files[0] ? evt.currentTarget.files[0].name : props.output && props.output.geneListFile ? props.output.geneListFile.url : ''
                setGeneListFile(({ url: _, ...file }) => ({ ...file, url }));
              }}
              type="file"
              accept=".csv,.txt"
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
              <button className="btn btn-lg">Gene List File</button>
              <span className="text-lg">{geneListFile.filename||geneListFile.url||'No file chosen'}</span>
            </div>
            <div className="bp4-input-group">
              <input
                type="text"
                name="description"
                className="bp4-input"
                placeholder={`File description`}
                onChange={evt => {
                  setGeneListFile(({ description: _, ...file }) => ({ ...file, description: evt.target.value }))
                }}
                value={geneListFile.description||''}
              />
            </div>
            <div className="inline-flex flex-row">
              <a
                href="/api/v1/components/service/ctd/example_geneSetCTD.csv"
                download="example_geneSetCTD.csv"
              >
                <Bp4Button
                    large
                    text="Example"
                    rightIcon="bring-data"
                />
              </a>
            </div>

            <p style={{ fontSize: "16px", marginTop: "20px" }}><b>Select Adj. Matrix file</b>, acceptable file type: <b>.csv</b></p>
            <input
              ref={adjMatrixInputRef}
              name="adjMatrixFile"
              className={classNames('file-input file-input-lg', { 'hidden': adjMatrixFile.url })}
              onChange={evt => {
                const url = evt.currentTarget && evt.currentTarget.files && evt.currentTarget.files[0] ? evt.currentTarget.files[0].name : props.output && props.output.adjMatrixFile ? props.output.adjMatrixFile.url : ''
                setAdjMatrixFile(({ url: _, ...file }) => ({ ...file, url }));
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
              <button className="btn btn-lg">Adj. Matrix File</button>
              <span className="text-lg">{adjMatrixFile.filename||adjMatrixFile.url||'No file chosen'}</span>
            </div>
            <div className="bp4-input-group">
              <input
                type="text"
                name="description"
                className="bp4-input"
                placeholder={`File description`}
                onChange={evt => {
                  setAdjMatrixFile(({ description: _, ...file }) => ({ ...file, description: evt.target.value }))
                }}
                value={adjMatrixFile.description||''}
              />
            </div>
            <div className="inline-flex flex-row">
              <a
                href="/api/v1/components/service/ctd/example_geneSetCTD.csv"
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
                  (props.output && props.output.geneListFile && (geneListFile.url === undefined || (geneListFile.description === props.output.geneListFile.description && geneListFile.url === props.output.geneListFile.url))) ||
                  (props.output && props.output.adjMatrixFile && (adjMatrixFile.url === undefined || (adjMatrixFile.description === props.output.adjMatrixFile.description && adjMatrixFile.url === props.output.adjMatrixFile.url)))
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
  .story(props => "Files uploaded")
  .build()

export const CTDUseCustomMatrixFileURLs = MetaNode('CTDUseCustomMatrixFileURLs')
  .meta({
    label: 'CTDUseCustomMatrixFileURLs',
    description: 'URLs to the files for CTD Cusotm Matrix Use API call',
    icon: [file_icon],
  })
  .codec(z.object({
    geneListFile: FileC,
    adjMatrixFile: FileC,
    rDataFile: FileC,
  }))
  .view(fileURLs => {
    const filenames = ['geneListFile', 'adjMatrixFile', 'rDataFile'] as const
    return (
      <Table
      height={500}
      cellRendererDependencies={[fileURLs]}
      numRows={filenames.length}
      enableGhostCells
      enableFocusedCell
      >
      <Column
        name="File Name"
        cellRenderer={row => <Cell key={row+''}>{fileURLs[filenames[row]].filename}</Cell>}
      />
      <Column
        name="URL"
        cellRenderer={row => <Cell key={row+''}>{fileURLs[filenames[row]].url}</Cell>}
      />
      <Column
        name="Size"
        cellRenderer={row => <Cell key={row+''}>{fileURLs[filenames[row]].size}</Cell>}
      />
      <Column
        name="Description"
        cellRenderer={row => <Cell key={row+''}>{fileURLs[filenames[row]].description}</Cell>}
      />
      </Table>
    )
  })
  .build()


export const CTDUseCustomMatrixFileInput = MetaNode('CTDUseCustomMatrixFileInput')
  .meta({
    label: 'CTD - UseCustomMatrix',
    description: 'Upload three files, Gene list, Adj. Matrix and RData, for CTD custom matrix.',
    icon: [input_icon],
  })
  .inputs()
  .output(CTDUseCustomMatrixFileURLs)
  .prompt(props => { 
    const { data: session } = useSessionWithId()

    const [geneListFile, setGeneListFile] = React.useState<z.infer<typeof FileC>>({ url: '', filename: '' })
    const geneListInputRef = React.useRef<HTMLInputElement>(null)

    const [adjMatrixFile, setAdjMatrixFile] = React.useState<z.infer<typeof FileC>>({ url: '', filename: '' })
    const adjMatrixInputRef = React.useRef<HTMLInputElement>(null)

    const [rDataFile, setRDataFileFile] = React.useState<z.infer<typeof FileC>>({ url: '', filename: '' })
    const rDataInputRef = React.useRef<HTMLInputElement>(null)

    React.useEffect(() => {
      if (!props.output) return
      setGeneListFile(props.output.geneListFile)
      setAdjMatrixFile(props.output.adjMatrixFile)
      setRDataFileFile(props.output.rDataFile)
    }, [props.output])
    return (
      <div>
        {!session || !session.user ? (
          <div className="alert alert-warning shadow-lg block prose">
            You are required to &nbsp; <button className="btn btn-sm" onClick={() => {Auth.signIn()}}>sign in</button> &nbsp; to upload files.
          </div>
        ) : (
          <form
            className="my-2 inline-flex flex-col"
            onSubmit={async (evt) => {
              evt.preventDefault();
              if (props.output != null && geneListFile.url && geneListFile.url === props.output.geneListFile?.url && geneListFile.filename && geneListFile.size &&
                        adjMatrixFile.url && adjMatrixFile.url === props.output.adjMatrixFile?.url && adjMatrixFile.filename && adjMatrixFile.size) {                
                props.submit({ geneListFile, adjMatrixFile, rDataFile })
              }else{
                const formData = new FormData(evt.currentTarget)
                formData.delete("description");
                let response = await clientUploadFile(formData);   
                const geneSetFileResponse = response.geneSetFile[0];
                const adjMatrixFileResponse = response.adjMatrixFile[0];   
                const rDataFileResponse = response.rDataFile[0];            

                props.submit({
                  geneListFile: {
                    description: "Gene list file",
                    url: geneSetFileResponse.url,
                    filename: geneSetFileResponse.filename,
                    size: geneSetFileResponse.size,
                  },
                  adjMatrixFile: {
                    description: "Adj. Matrix File",
                    url: adjMatrixFileResponse.url,
                    filename: adjMatrixFileResponse.filename,
                    size: adjMatrixFileResponse.size,
                  },
                  rDataFile: {
                    description: "RData File",
                    url: rDataFileResponse.url,
                    filename: rDataFileResponse.filename,
                    size: rDataFileResponse.size,
                  }
                })
              }
            }}
          >
            <p style={{ fontSize: "16px" }}><b>Select gene list file</b>, acceptable file types: <b>.csv, .txt</b></p>
            <input
              ref={geneListInputRef}
              name="geneSetFile"
              className={classNames('file-input file-input-lg', { 'hidden': geneListFile.url })}
              onChange={evt => {
                const url = evt.currentTarget && evt.currentTarget.files && evt.currentTarget.files[0] ? evt.currentTarget.files[0].name : props.output && props.output.geneListFile ? props.output.geneListFile.url : ''
                setGeneListFile(({ url: _, ...file }) => ({ ...file, url }));
              }}
              type="file"
              accept=".csv,.txt"
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
              <button className="btn btn-lg">Gene List File</button>
              <span className="text-lg">{geneListFile.filename||geneListFile.url||'No file chosen'}</span>
            </div>
            <div className="bp4-input-group">
              <input
                type="text"
                name="description"
                className="bp4-input"
                placeholder={`File description`}
                onChange={evt => {
                  setGeneListFile(({ description: _, ...file }) => ({ ...file, description: evt.target.value }))
                }}
                value={geneListFile.description||''}
              />
            </div>
            <div className="inline-flex flex-row">
              <a
                href="/api/v1/components/service/ctd/example_geneSetCTD.csv"
                download="example_geneSetCTD.csv"
              >
                <Bp4Button
                    large
                    text="Example"
                    rightIcon="bring-data"
                />
              </a>
            </div>

            <p style={{ fontSize: "16px", marginTop: "20px" }}><b>Select Adj. Matrix file</b>, acceptable file type: <b>.csv</b></p>
            <input
              ref={adjMatrixInputRef}
              name="adjMatrixFile"
              className={classNames('file-input file-input-lg', { 'hidden': adjMatrixFile.url })}
              onChange={evt => {
                const url = evt.currentTarget && evt.currentTarget.files && evt.currentTarget.files[0] ? evt.currentTarget.files[0].name : props.output && props.output.adjMatrixFile ? props.output.adjMatrixFile.url : ''
                setAdjMatrixFile(({ url: _, ...file }) => ({ ...file, url }));
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
              <button className="btn btn-lg">Adj. Matrix File</button>
              <span className="text-lg">{adjMatrixFile.filename||adjMatrixFile.url||'No file chosen'}</span>
            </div>
            <div className="bp4-input-group">
              <input
                type="text"
                name="description"
                className="bp4-input"
                placeholder={`File description`}
                onChange={evt => {
                  setAdjMatrixFile(({ description: _, ...file }) => ({ ...file, description: evt.target.value }))
                }}
                value={adjMatrixFile.description||''}
              />
            </div>
            <div className="inline-flex flex-row">
              <a
                href="/api/v1/components/service/ctd/example_geneSetCTD.csv"
                download="example_geneSetCTD.csv"
              >
                <Bp4Button
                    large
                    text="Example"
                    rightIcon="bring-data"
                />
              </a>
            </div>


            <p style={{ fontSize: "16px", marginTop: "20px" }}><b>Select custom RData file</b>, acceptable file type: <b>.RData</b></p>
            <input
              ref={rDataInputRef}
              name="rDataFile"
              className={classNames('file-input file-input-lg', { 'hidden': rDataFile.url })}
              onChange={evt => {
                const url = evt.currentTarget && evt.currentTarget.files && evt.currentTarget.files[0] ? evt.currentTarget.files[0].name : props.output && props.output.rDataFile ? props.output.rDataFile.url : ''
                setRDataFileFile(({ url: _, ...file }) => ({ ...file, url }));
              }}
              type="file"
              accept=".RData"
            />
            <div
              className={classNames('inline-flex flex-row items-center gap-4', { 'hidden': !rDataFile.url })}
              onClick={evt => {
                evt.preventDefault()
                if (rDataInputRef.current) {
                  rDataInputRef.current.click()
                }
              }}
            >
              <button className="btn btn-lg">Adj. Matrix File</button>
              <span className="text-lg">{rDataFile.filename||rDataFile.url||'No file chosen'}</span>
            </div>
            <div className="bp4-input-group">
              <input
                type="text"
                name="description"
                className="bp4-input"
                placeholder={`File description`}
                onChange={evt => {
                  setRDataFileFile(({ description: _, ...file }) => ({ ...file, description: evt.target.value }))
                }}
                value={rDataFile.description||''}
              />
            </div>
            <div className="inline-flex flex-row">
              <a
                href="/api/v1/components/service/ctd/example_geneSetCTD.csv"
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
                  (props.output && props.output.geneListFile && (geneListFile.url === undefined || (geneListFile.description === props.output.geneListFile.description && geneListFile.url === props.output.geneListFile.url))) ||
                  (props.output && props.output.adjMatrixFile && (adjMatrixFile.url === undefined || (adjMatrixFile.description === props.output.adjMatrixFile.description && adjMatrixFile.url === props.output.adjMatrixFile.url))) ||
                  (props.output && props.output.rDataFile && (rDataFile.url === undefined || (rDataFile.description === props.output.rDataFile.description && rDataFile.url === props.output.rDataFile.url)))
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
  .story(props => "Files uploaded")
  .build()
