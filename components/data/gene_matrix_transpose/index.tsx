import { MetaNode } from '@/spec/metanode'
import { z } from 'zod'
import * as dict from '@/utils/dict'
import * as array from '@/utils/array'
import { Table2 as Table, Cell, Column } from '@blueprintjs/table'
import { GeneSet } from '@/components/core/input/set'

export const GMT = MetaNode.createData(`GMT`)
  .meta({
    label: `Gene Matrix Transpose`,
    description: 'Terms mapped to genes',
    example: {'A': ['B', 'C', 'D'], 'F': ['G', 'H']},
  })
  .codec(z.record(z.string(), z.array(z.string())))
  .view(gmt => {
    const gmt_items = dict.items(gmt)
    return (
      <div style={{ height: 500 }}>
        <Table
          cellRendererDependencies={[gmt_items]}
          numRows={gmt_items.length}
          enableGhostCells
          enableFocusedCell
        >
          <Column
            name="Term"
            cellRenderer={row => <Cell key={row+''}>{gmt_items[row].key.toString()}</Cell>}
          />
          <Column
            name="Geneset"
            cellRenderer={row => <Cell key={row+''}>{gmt_items[row].value.join('\t')}</Cell>}
          />
        </Table>
      </div>
    )
  })
  .build()

export const GMTUnion = MetaNode.createProcess('GMTUnion')
  .meta({
    label: `Compute Union Geneset`,
    description: 'Find the union set of all genes in the GMT'
  })
  .inputs({ gmt: GMT })
  .output(GeneSet)
  .resolve(async (props) => {
    return array.unique(dict.values(props.inputs.gmt).flatMap(geneset => geneset))
  })
  .build()

export const GMTIntersection = MetaNode.createProcess('GMTIntersection')
  .meta({
    label: `Compute Intersection Geneset`,
    description: 'Find the intersecting set of all genes in the GMT'
  })
  .inputs({ gmt: GMT })
  .output(GeneSet)
  .resolve(async (props) => {
    return dict.values(props.inputs.gmt).reduce((A, B) => array.intersection(A, B))
  })
  .build()

export const GMTConsensus = MetaNode.createProcess('GMTConsensus')
  .meta({
    label: `Compute Consensus Geneset`,
    description: 'Find genes which appear in more than one set'
  })
  .inputs({ gmt: GMT })
  .output(GeneSet)
  .resolve(async (props) => {
    const gene_counts: Record<string, number> = {}
    dict.values(props.inputs.gmt)
      .forEach(geneset =>
        geneset.forEach(gene =>
          gene_counts[gene] = (gene_counts[gene]||0)+1
        )
      )
    return dict.items(gene_counts)
      .filter(({ value }) => value > 1)
      .map(({ key }) => key as string)
  })
  .build()
