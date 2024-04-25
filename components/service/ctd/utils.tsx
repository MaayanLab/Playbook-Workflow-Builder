import { MetaNode } from '@/spec/metanode'
import { z } from 'zod'
import { FilePrompt, FileC } from '@/components/core/file'
import { downloadUrl } from '@/utils/download'
import { Table, Cell, Column} from '@/app/components/Table'
import SafeRender from '@/utils/saferender'
import {datafile_icon, ctd_icon, file_icon } from '@/icons'
import { clientLoadExample } from  '@/components/data/anndata/api/example.h5ad/client'

export const RData = MetaNode('RData')
.meta({
  label: 'CTD RData',
  description: 'An RData file',
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
    <div>
      <p><b>RData</b></p>
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
    </div>
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
    <div>
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
    </div>
  )
})
.build()

export const AdjacencyMatrixFileUpload = MetaNode('AdjacencyMatrixFileUpload')
.meta({
  label: 'Upload a CTD Adjacency Matrix',
  description: 'A file containing an Adjacency Matrix for CTD - "Connect the Dots in Precalculated Graph"',
  tags: {
    Type: {
      File: 1,
      CTD: 1,
      Matrix: 1
    }
  },
  icon: [file_icon]
})
.codec(FileC)
.inputs()
.output(AdjacencyMatrix)
.prompt(props => <><FilePrompt {...props} example={clientLoadExample} />{props.output ? <SafeRender component={AdjacencyMatrix.view} props={props.output} /> : null}</>)
.resolve(async (props) => {
  return {
    url: props.data.url,
    filename: props.data.filename
  };
}).story(props =>
  `An Adjacency matrix${props.data && props.data.description ? ` containing ${props.data.description}` : ''} was uploaded.`
)
.build()

export const GeneExpressions = MetaNode('GeneExpressions')
.meta({
  label: 'Gene Expressions',
  description: 'A Gene Expressions file',
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
    <div>
    <p><b>Gene Expressions</b></p>
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
    </div> 
  )
})
.build()

export const GeneExpressionsFileUpload = MetaNode('GeneExpressionsFileUpload')
.meta({
  label: 'Upload a Gene Expressions File',
  description: 'A file containing a list of Genes that can generate a CTD Adjacency Matrix for a custom CTD Precalculations ("Connect the Dots in Precalculated Graph")',
  tags: {
    Type: {
      File: 1,
      Gene: 1,
      CTD: 1
    }
  },
  icon: [file_icon]
})
.codec(FileC)
.inputs()
.output(GeneExpressions)
.prompt(props => <><FilePrompt {...props} example={clientLoadExample} />{props.output ? <SafeRender component={GeneExpressions.view} props={props.output} /> : null}</>)
.resolve(async (props) => {
  return {
    url: props.data.url,
    filename: props.data.filename
  };
}).story(props =>
  `A  Gene Expressions file ${props.data && props.data.description ? ` containing ${props.data.description}` : ''} was uploaded.`
)
.build()

export const CTD_FileDownload = MetaNode('CTD_FileDownload')
  .meta({
    label: 'CTD File Download',
    description: '',
    icon: [datafile_icon]
  })
  .codec(
    z.object({
      url: z.string(),
      filename: z.string(),
    })
  )
  .view(props => {
    return (
      <div>
        <p><b>CTD File Download</b></p>
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
      </div>
    )
  }).build()
