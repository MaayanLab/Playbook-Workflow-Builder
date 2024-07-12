import React from 'react'
import { MetaNode } from "@/spec/metanode"
import { z } from 'zod'
import { AnnData } from '@/components/data/anndata'
import { useClientMetadataFromFile } from './api/metadata/client'
import { updateMetadataColumn } from './api/metadata/update'
import classNames from 'classnames'
import { GeneCountMatrix } from '../gene_count_matrix'

export const LabelAnnDataMetadata = MetaNode('LabelAnnDataMetadata')
  .meta({
    label: 'Manual Labeling of Samples',
    description: 'Manually label samples as either control or perturbation.',
  })
  .codec(z.record(z.string(), z.record(z.string(), z.string())))
  .inputs({ matrix: AnnData })
  .output(AnnData)
  .prompt(props => {
    const matrix = props.output ? props.output : props.inputs.matrix
    const { data, isLoading } = useClientMetadataFromFile(matrix)
    const [tableData, setTableData] = React.useState(props.data ? props.data : {} as Record<string, Record<string, string>>)

    React.useEffect(() => {
      if (!data) return
      setTableData(data)
    }, [data])

    const { index, columns } = React.useMemo(() => {
      const index: Record<string, undefined> = {}
      const columns: string[] = []
      if (tableData) {
        Object.keys(tableData).forEach(column => {
          columns.push(column)
          Object.keys(tableData[column]).forEach(ind => {
            index[ind] = undefined
          })
        })
      }
      if (!columns.includes('Type: Control or Perturbation')) {
        columns.push('Type: Control or Perturbation')
      }
      return { index: Object.keys(index), columns }
    }, [tableData])

    const handleCellChange = React.useCallback((row: string, column: string, value: string) => {
      const newTableData = {...tableData}
      if (!(column in newTableData)) newTableData[column] = {}
      newTableData[column][row] = value
      setTableData(newTableData)
    }, [tableData])

    return (
      <div>
        <progress className={classNames('progress w-full', { 'hidden': !isLoading })}></progress>
        <div className={classNames("overflow-x-auto", { 'hidden': index.length === 0 || columns.length === 0 })}>
          <table style={{ borderCollapse: 'collapse', border: '1px solid black' }}>
            <thead>
              <tr>
                <th style={{ border: '1px solid black', padding: '5px' }}>&nbsp;</th>
                {columns.map((column, columnIndex) => (
                  <th key={columnIndex} style={{ border: '1px solid black', padding: '5px' }}>{column}</th>
                ))}
              </tr>
            </thead>
            <tbody>
              {index.map((row, rowIndex) => (
                <tr key={rowIndex}>
                  <td
                    style={{
                      border: '1px solid black',
                      padding: '5px',
                    }}
                  >{row}</td>
                  {columns.map((column, columnIndex) => (
                    <td
                      key={columnIndex}
                      style={{
                        border: '1px solid black',
                        padding: '5px',
                      }}
                    >
                      {column === 'Type: Control or Perturbation' ?
                      <div className="flex justify-center place-items-center gap-2">
                        <span className="prose max-w-none hover:cursor-pointer" onClick={() => {
                          handleCellChange(row, column, 'Control')
                        }}>Control</span> 
                        <input
                          type="checkbox"
                          className="toggle"
                          ref={ref => {if (ref) ref.indeterminate = !['Control', 'Perturbation'].includes((tableData[column] || {})[row])}}
                          checked={(tableData[column] || {})[row] === 'Perturbation'}
                          onClick={evt =>
                            handleCellChange(row, column, (tableData[column] || {})[row] === 'Control' ? 'Perturbation' : 'Control')
                          }
                        />
                        <span className="prose max-w-none hover:cursor-pointer" onClick={() => {
                          handleCellChange(row, column, 'Perturbation')
                        }}>Perturbation</span>
                      </div>
                      : <input
                        type="text"
                        value={tableData[column][row]}
                        onChange={evt =>
                          handleCellChange(row, column, evt.target.value)
                        }
                        style={{
                          width: '100%',
                          boxSizing: 'border-box',
                          border: 'none',
                          outline: 'none',
                          textAlign: 'center',
                        }}
                        />}
                    </td>
                  ))}
                </tr>
              ))}
            </tbody>
          </table>
        </div>
        <button className="bp5-button bp5-large" onClick={() => {
          props.submit(tableData, true)
        }}>Submit</button>
      </div>
    );
  })
  .resolve(async (props) => {
    return await updateMetadataColumn({
      file: props.inputs.matrix,
      data: props.data,
    })
  })
  .story(props => ({
    abstract: 'The samples were then labeled as either control or perturbation to allow for further analysis.',
  }))
  .build()

export const LabelGeneCountMatrix = MetaNode('LabelGeneCountMatrix')
  .meta(LabelAnnDataMetadata.meta)
  .codec(LabelAnnDataMetadata.codec)
  .inputs({ matrix: GeneCountMatrix })
  .output(LabelAnnDataMetadata.output)
  .prompt(LabelAnnDataMetadata.prompt)
  .resolve(LabelAnnDataMetadata.resolve)
  .story(LabelAnnDataMetadata.story)
  .build()
