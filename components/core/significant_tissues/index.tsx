import React from 'react'
import { MetaNode } from '@/spec/metanode'
import { z } from 'zod'
import { tissue_icon, weighted_icon } from '@/icons'
import { Table2 as Table, Column, Cell } from '@blueprintjs/table'

export const SignificantTissueC = z.array(z.object({
  tissue: z.string().describe('The tissue term'),
  zscore: z.number().describe('A zscore for the significance'),
}))

export const SignificantTissues = MetaNode.createData('SignificantTissues')
  .meta({
    label: 'SignificantTissues',
    description: 'Tissues scored using a combined stouffer statistic',
    icon: [weighted_icon, tissue_icon],
  })
  .codec(SignificantTissueC)
  .view(tissues => {
    return (
      <div style={{ height: 500 }}>
        <Table
          cellRendererDependencies={[tissues]}
          numRows={tissues.length}
          enableGhostCells
          enableFocusedCell
        >
          <Column
            name="Tissue"
            cellRenderer={row => <Cell key={row+''}>{tissues[row].tissue}</Cell>}
          />
          <Column
            name="ZScore"
            cellRenderer={row => <Cell key={row+''}>{tissues[row].zscore.toPrecision(3)}</Cell>}
          />
        </Table>
      </div>
    )
  })
  .build()
