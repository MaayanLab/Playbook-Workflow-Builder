import React from 'react'
import { MetaNode } from '@/spec/metanode'
import { Primative, TumorGeneExpression} from '@/components/core/input/primitives'
import { z } from 'zod'
import { Table, Cell, Column } from '@/app/components/Table'
import * as array from '@/utils/array'



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
  })
  .codec(z.array(z.object({ Disease : z.string(), 
                            Gene_symbol : z.string(),
                            TPM_mean : z.number()})))
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
          name={"Disease"}
          cellRenderer={row => <Cell key={row+''}>{expressionTable[row].Disease}</Cell>}
        />
        <Column
          name="GeneSymbol"
          cellRenderer={row => <Cell key={row+''}>{expressionTable[row].Gene_symbol}</Cell>}
        />
        <Column
          name="TpmMean"
          cellRenderer={row => <Cell key={row+''}>{expressionTable[row].TPM_mean}</Cell>}
        />
      </Table>
    )
  })
  .build()


export const GeneExpressionInTumor = GeneExpressionInTumor_T(TumorGeneExpression)