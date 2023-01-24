import React from 'react'
import { MetaNode } from '@/spec/metanode'
import { z } from 'zod'
import { Gene, Drug, Primative, Pathway, Phenotype, Tissue, Disease } from '@/components/core/input/primitives'
import { Table2 as Table, Column, Cell } from '@blueprintjs/table'
import { weighted_icon } from '@/icons'
import * as array from '@/utils/array'

const Scored_T = (T: Primative) => MetaNode.createData(`Scored[${T.name}]`)
  .meta({
    label: `Scored ${T.label}s`,
    description: `ZScores of ${T.label}s`,
    color: T.color,
    icon: [...array.ensureArray(T.icon), weighted_icon],
    tags: {
      Type: {
        [T.label]: 1,
      },
      Cardinality: {
        Scored: 1,
      },
    },
  })
  .codec(z.array(z.object({ term: z.string(), zscore: z.number() })))
  .view(scored => {
    return (
      <div style={{ height: 500 }}>
        <Table
          cellRendererDependencies={[scored]}
          numRows={scored.length}
          enableGhostCells
          enableFocusedCell
        >
          <Column
            name={T.label}
            cellRenderer={row => <Cell key={row+''}>{scored[row].term}</Cell>}
          />
          <Column
            name="ZScore"
            cellRenderer={row => <Cell key={row+''}>{scored[row].zscore.toPrecision(3)}</Cell>}
          />
        </Table>
      </div>
    )
  })
  .build()

export const ScoredDiseases = Scored_T(Disease)
export const ScoredDrugs = Scored_T(Drug)
export const ScoredGenes = Scored_T(Gene)
export const ScoredPathways = Scored_T(Pathway)
export const ScoredPhenotypes = Scored_T(Phenotype)
export const ScoredTissues = Scored_T(Tissue)
