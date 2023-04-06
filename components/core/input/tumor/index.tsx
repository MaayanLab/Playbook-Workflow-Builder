import React from 'react'
import { MetaNode } from '@/spec/metanode'
import { Primative, TumorGeneExpression} from '@/components/core/input/primitives'
import { z } from 'zod'
import { Table, Cell, Column } from '@/app/components/Table'
import * as array from '@/utils/array'
import { tumor_icon } from '@/icons'



const GeneExpressionInTumor_T = (T: Primative) => MetaNode(`[${T.name}]`)
  .meta({
    label: `${T.label}`,
    description: `Output of ${T.label}s`,
    color: T.color,
    tags: {
      Type: {
        [T.label]: 1,
      },
      Cardinality: {
        Scored: 1,
      },
    },
    icon: [tumor_icon],
  })
  .codec(z.array(z.object({ TPM_mean : z.number(), 
                            TPM_sd : z.number(), 
                            TPM_median : z.number(),
                            Disease : z.string(), 
                            Gene_symbol : z.string(),
                            Gene_Ensembl_ID : z.string(),
                            Dataset : z.string()})))
  .view(expressionTable => {
    return (
      <Table
        height={500}
        cellRendererDependencies={[expressionTable]}
        numRows={expressionTable.length}
        enableGhostCells
        enableFocusedCell
      > 
        <Column
          name="Gene ID"
          cellRenderer={row => <Cell key={row+''}>{expressionTable[row].Gene_symbol}</Cell>}
        />
        <Column
          name="Data Set"
          cellRenderer={row => <Cell key={row+''}>{expressionTable[row].Dataset}</Cell>}
        />
        <Column
          name={"Disease"}
          cellRenderer={row => <Cell key={row+''}>{expressionTable[row].Disease}</Cell>}
        />
        <Column
          name={"TPM Mean"}
          cellRenderer={row => <Cell key={row+''}>{expressionTable[row].TPM_mean}</Cell>}
        />
        <Column
          name={"TPM Stand. Dev."}
          cellRenderer={row => <Cell key={row+''}>{expressionTable[row].TPM_sd}</Cell>}
        />
        <Column
          name={"TPM Median"}
          cellRenderer={row => <Cell key={row+''}>{expressionTable[row].TPM_median}</Cell>}
        />
        <Column
          name="Ensembl ID"
          cellRenderer={row => <Cell key={row+''}>{expressionTable[row].Gene_Ensembl_ID}</Cell>}
        />
      </Table>
    )
  })
  .build()


export const GeneExpressionInTumor = GeneExpressionInTumor_T(TumorGeneExpression)