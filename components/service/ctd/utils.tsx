import { MetaNode } from '@/spec/metanode'
import { z } from 'zod'
import { downloadUrl } from '@/utils/download'
import { Table, Cell, Column} from '@/app/components/Table'
import {datafile_icon } from '@/icons'


export const CTD_DataSetC = z.object({
  ctdPermutations:z.object({
    url: z.string(),
    filename: z.string(),
  }),
  ctdMatrix: z.object({
    url: z.string(),
    filename: z.string(),
  }),
  geneSet: z.array(z.string())
});


export const CTD_DataSet = MetaNode('CTD_DataSet')
.meta({
  label: 'Data set containg all of the inout types for the final CTD request!',
  description: 'Data set containig the Adjacency Matrix (.csv), PErmutations (.RData) and Gene Set (.csv)',
  icon: [datafile_icon],
})
.codec(CTD_DataSetC)
.view(props => {
  return (
    <>
      <p><b>CTD Permutations</b></p>
      <Table
      cellRendererDependencies={[1]}
      numRows={1}
      downloads={{
        'URL': () => downloadUrl(props.ctdPermutations.url, props.ctdPermutations.filename)
      }}
      >
      <Column
        name="Filename"
        cellRenderer={row => <Cell key={row+''}>{props.ctdPermutations.filename}</Cell>}
      />
      <Column
        name="URL"
        cellRenderer={row => <Cell key={row+''}>{props.ctdPermutations.url}</Cell>}
      />
      </Table>

      <p><b>CTD Matrix</b> (You can dowload this file in the "CTD Adjacency Matrix" output card)</p>
      <Table
      cellRendererDependencies={[1]}
      numRows={1}
      >
      <Column
        name="Filename"
        cellRenderer={row => <Cell key={row+''}>{props.ctdMatrix.filename}</Cell>}
      />
      <Column
        name="URL"
        cellRenderer={row => <Cell key={row+''}>{props.ctdMatrix.url}</Cell>}
      />
      </Table>

      <p><b>Gene Set</b></p>
      <div style={{ overflowY: 'auto', height: '300px',  width: '300px' }}>{( props.geneSet).join('\n')}</div>
    </>
  )
})
.build()

export const AdjacencyMatrix = MetaNode('AdjacencyMatrix')
.meta({
  label: 'CTD Adjacency Matrix',
  description: 'An Adjacency matrix file',
  icon: [datafile_icon],
})
.codec(
  z.object({
    url: z.string(),
    filename: z.string(),
  })
)
.view(props => {
  return (
    <>
      <p><b>Adjacency Matrix</b></p>
      <Table
      cellRendererDependencies={[1]}
      numRows={1}
      downloads={{
        'URL': () => downloadUrl(props.url, props.filename)
      }}
      >
      <Column
        name="Filename"
        cellRenderer={row => <Cell key={row+''}>{props.filename}</Cell>}
      />
      <Column
        name="URL"
        cellRenderer={row => <Cell key={row+''}>{props.url}</Cell>}
      />
      </Table>
    </>
  )
})
.build()

