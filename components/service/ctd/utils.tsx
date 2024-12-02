import { MetaNode } from '@/spec/metanode'
import { z } from 'zod'
import { downloadUrl } from '@/utils/download'
import { Table, Cell, Column} from '@/app/components/Table'
import {datafile_icon } from '@/icons'


export const CTD_DataSetC = z.object({
  ctdPermutations:z.object({
    url: z.string(),
    filename: z.string(),
    size: z.number()
  }),
  expressions: z.object({
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
  external: true,
})
.codec(CTD_DataSetC)
.view(props => {
  return (
    <>
      <p><b>CTD Permutations</b> (JSON)</p>
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
      <Column
        name="Size"
        cellRenderer={row => <Cell key={row+''}>{ formatBytes(props.ctdPermutations.size) }</Cell>}
      />
      </Table>

      <p><b>Gene Expressions</b> (CSV)</p>
      <Table
      cellRendererDependencies={[1]}
      numRows={1}
      >
      <Column
        name="Filename"
        cellRenderer={row => <Cell key={row+''}>{props.expressions.filename}</Cell>}
      />
      <Column
        name="URL"
        cellRenderer={row => <Cell key={row+''}>{props.expressions.url}</Cell>}
      />
      </Table>

      <p><b>Gene Set</b></p>
      <div style={{ overflowY: 'auto', height: '300px',  width: '300px' }}>{( props.geneSet).join('\n')}</div>
    </>
  )
})
.build()

export const CTDAdjacencyAndExpressions = MetaNode('CTDAdjacencyAndExpressions')
.meta({
  label: 'CTD Adjacency and Expression file',
  description: 'An Adjacency JSON file and Expression CSV file.',
  icon: [datafile_icon],
  external: true,
})
.codec(
  z.object({
    adjacency: z.object({
      url: z.string(),
      filename: z.string(),
      size: z.number()
    }),
    expressions: z.object({
      url: z.string(),
      filename: z.string()
    })
  })
)
.view(props => {
  return (
    <>
      <p><b>Adjacency JSON file</b></p>
      <Table
      cellRendererDependencies={[1]}
      numRows={1}
      downloads={{
        'URL': () => downloadUrl(props.adjacency.url, props.adjacency.filename)
      }}
      >
      <Column
        name="Filename"
        cellRenderer={row => <Cell key={row+''}>{props.adjacency.filename}</Cell>}
      />
      <Column
        name="URL"
        cellRenderer={row => <Cell key={row+''}>{props.adjacency.url}</Cell>}
      />
      <Column
        name="Size"
        cellRenderer={row => <Cell key={row+''}>{ formatBytes(props.adjacency.size) }</Cell>}
      />
      </Table>

      <p><b>Expressions CSV file</b></p>
      <Table
      cellRendererDependencies={[1]}
      numRows={1}
      downloads={{
        'URL': () => downloadUrl(props.expressions.url, props.expressions.filename)
      }}
      >
      <Column
        name="Filename"
        cellRenderer={row => <Cell key={row+''}>{props.expressions.filename}</Cell>}
      />
      <Column
        name="URL"
        cellRenderer={row => <Cell key={row+''}>{props.expressions.url}</Cell>}
      />
      </Table>
    </>
  )
})
.build()

function formatBytes(bytes: number, decimals = 2) {
  if (!+bytes) return '0 Bytes';

  const k = 1024;
  const dm = decimals < 0 ? 0 : decimals;
  const sizes = ['Bytes', 'KiB', 'MiB', 'GiB', 'TiB', 'PiB', 'EiB', 'ZiB', 'YiB'];

  const i = Math.floor(Math.log(bytes) / Math.log(k));

  return `${parseFloat((bytes / Math.pow(k, i)).toFixed(dm))} ${sizes[i]}`;
}

