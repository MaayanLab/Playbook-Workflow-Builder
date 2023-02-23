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
  .codec(z.array(z.object({ entrezGeneId : z.number(), 
                            patientId: z.string(),
                            sampleId: z.string(),
                            mutationType: z.string(),
                            mutationStatus: z.string(),
                            startPosition: z.number(),
                            endPosition: z.number(),
                            cancerTypeId: z.string()})))
  .view(scored => {
    return (
      <Table
        height={500}
        cellRendererDependencies={[scored]}
        numRows={scored.length}
        enableGhostCells
        enableFocusedCell
      >
        <Column
          name={"GeneId"}
          cellRenderer={row => <Cell key={row+''}>{scored[row].entrezGeneId}</Cell>}
        />
        <Column
          name="PatientId"
          cellRenderer={row => <Cell key={row+''}>{scored[row].patientId}</Cell>}
        />
        <Column
          name="SampleId"
          cellRenderer={row => <Cell key={row+''}>{scored[row].sampleId}</Cell>}
        />
        <Column
          name={"MutationType"}
          cellRenderer={row => <Cell key={row+''}>{scored[row].mutationType}</Cell>}
        />
        <Column
          name="MutationStatus"
          cellRenderer={row => <Cell key={row+''}>{scored[row].mutationStatus}</Cell>}
        />
        <Column
          name="StartPosition"
          cellRenderer={row => <Cell key={row+''}>{scored[row].startPosition}</Cell>}
        />
        <Column
          name="EndPosition"
          cellRenderer={row => <Cell key={row+''}>{scored[row].endPosition}</Cell>}
        />
        <Column
          name="CancerType"
          cellRenderer={row => <Cell key={row+''}>{scored[row].cancerTypeId}</Cell>}
        />
      </Table>
    )
  })
  .build()


export const GeneExpressionInTumor = GeneExpressionInTumor_T(TumorGeneExpression)