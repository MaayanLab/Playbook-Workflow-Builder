import React from 'react'
import { MetaNode } from '@/spec/metanode'
import { z } from 'zod'
import {
  Disease,
  Drug,
  Gene,
  Metabolite,
  Pathway,
  Phenotype,
  Primative,
  Protein,
  Glycan,
  RegulatoryElement,
  Tissue,
  Variant,
} from '@/components/core/primitives'
import { Table, Cell, Column } from '@/app/components/Table'
import { list_icon } from '@/icons'
import * as array from '@/utils/array'
import { downloadBlob } from '@/utils/download'
import pluralize from 'pluralize'
import {
  DiseaseSet,
  DrugSet,
  GeneSet,
  GlycanSet,
  MetaboliteSet,
  PathwaySet,
  PhenotypeSet,
  ProteinSet,
  RegulatoryElementSet,
  TissueSet,
  VariantSet,
} from '@/components/core/set'

const Ranked_T = (T: Primative) => MetaNode(`Ranked[${T.name}]`)
  .meta({
    label: `Ranked ${pluralize(T.label)}`,
    description: `From most significant to least`,
    icon: [...array.ensureArray(T.icon), list_icon],
    color: T.color,
    tags: {
      Type: {
        [T.label]: 1,
      },
      Cardinality: {
        Ranked: 1,
      },
    },
    ...(T.extra?.ranked?.meta || {}),
  })
  .codec(z.object({ description: z.string().optional(), ranked: z.array(z.string()) }))
  .view(({ ranked }) => {
    return (
      <Table
        height={500}
        cellRendererDependencies={[ranked]}
        numRows={ranked.length}
        enableGhostCells
        enableFocusedCell
        downloads={{
          JSON: () => downloadBlob(new Blob([JSON.stringify(ranked)], { type: 'application/json;charset=utf-8' }), 'data.json'),
          GMT: () => downloadBlob(new Blob([[`${T.label} Ranked`, '', ...ranked].join('\t')], { type: 'text/tab-separated-values;charset=utf-8' }), 'data.gmt'),
        }}
      >
        <Column
          name={T.label}
          cellRenderer={row => <Cell key={row+''}>{ranked[row]}</Cell>}
        />
      </Table>
    )
  })
  .build()

export const DiseaseRanked = Ranked_T(Disease)
export const DrugRanked = Ranked_T(Drug)
export const GeneRanked = Ranked_T(Gene)
export const ProteinRanked = Ranked_T(Protein)
export const GlycanRanked = Ranked_T(Glycan)
export const VariantRanked = Ranked_T(Variant)
export const RegulatoryElementRanked = Ranked_T(RegulatoryElement)
export const PathwayRanked = Ranked_T(Pathway)
export const PhenotypeRanked = Ranked_T(Phenotype)
export const TissueRanked = Ranked_T(Tissue)
export const MetaboliteRanked = Ranked_T(Metabolite)

export const RankedToSetT = [
  { RankedT: DiseaseRanked, SetT: DiseaseSet, },
  { RankedT: DrugRanked, SetT: DrugSet, },
  { RankedT: GeneRanked, SetT: GeneSet, },
  { RankedT: ProteinRanked, SetT: ProteinSet, },
  { RankedT: GlycanRanked, SetT: GlycanSet, },
  { RankedT: VariantRanked, SetT: VariantSet, },
  { RankedT: RegulatoryElementRanked, SetT: RegulatoryElementSet, },
  { RankedT: PathwayRanked, SetT: PathwaySet, },
  { RankedT: PhenotypeRanked, SetT: PhenotypeSet, },
  { RankedT: TissueRanked, SetT: TissueSet, },
  { RankedT: MetaboliteRanked, SetT: MetaboliteSet, },
].map(({ RankedT, SetT }) =>
  MetaNode(`RankedToSet[${RankedT.spec}]`)
    .meta({
      label: `Create Set from ${pluralize(RankedT.meta.label)}`,
      description: 'Ignore ordering of list'
    })
    .inputs({ ranked: RankedT })
    .output(SetT)
    .resolve(async (props) => {
      return { set: props.inputs.ranked.ranked, description: props.inputs.ranked.description }
    })
    .story(props => ``)
    .build()
)
