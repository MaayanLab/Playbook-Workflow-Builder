import React from 'react'
import { MetaNode } from '@/spec/metanode'
import { z } from 'zod'
import { supervenn_icon } from '@/icons'
import dynamic from 'next/dynamic'
import * as dict from '@/utils/dict'
import { GMT } from '@/components/data/gene_matrix_transpose'
import { DMT } from '@/components/data/drug_matrix_transpose'

const ReactSupervenn = dynamic(() => import('react-supervenn'), { ssr: false })

const hyposet_url = 'https://hyposet.maayanlab.cloud'

export const Supervenn = MetaNode('Supervenn')
  .meta({
    label: 'Supervenn Visualization',
    description: 'A visualization for comparing sets',
    icon: [supervenn_icon],
  })
  .codec(z.any())
  .view(props => (
    <div className="flex" style={{ minHeight: 500 }}>
      <ReactSupervenn {...props} />
    </div>
  ))
  .build()

export const SupervennFromGMT = MetaNode('SupervennFromGMT')
  .meta({
    label: 'Compare sets with Supervenn',
    description: 'Interactively analyse overlap between sets',
    icon: [supervenn_icon],
  })
  .inputs({ gmt: GMT })
  .output(Supervenn)
  .resolve(async (props) => {
    const inputProps = {
      sets: dict.values(props.inputs.gmt).map(({ set }) => set),
      set_annotations: dict.keys(props.inputs.gmt),
      widths_minmax_ratio: 0.1,
      rotate_col_annotations: true,
    }
    const res = await fetch(`${hyposet_url}/api/supervenn`, {
      method: 'POST',
      headers: {
        'Content-Type': 'application/json',
        'Accept': 'application/json',
      },
      body: JSON.stringify(inputProps),
    })
    return await res.json()
  })
  .story(props =>
    `The collection of gene sets was then visualized with a Supervenn diagram${''/* [FIGURE]*/}.`
  )
  .build()

export const SupervennFromDMT = MetaNode('SupervennFromDMT')
  .meta({
    label: 'Compare sets with Supervenn',
    description: 'Interactively analyse overlap between sets',
    icon: [supervenn_icon],
  })
  .inputs({ dmt: DMT })
  .output(Supervenn)
  .resolve(async (props) => {
    const inputProps = {
      sets: dict.values(props.inputs.dmt).map(({ set }) => set),
      set_annotations: dict.keys(props.inputs.dmt),
      widths_minmax_ratio: 0.1,
      rotate_col_annotations: true,
    }
    const res = await fetch(`${hyposet_url}/api/supervenn`, {
      method: 'POST',
      headers: {
        'Content-Type': 'application/json',
        'Accept': 'application/json',
      },
      body: JSON.stringify(inputProps),
    })
    return await res.json()
  })
  .story(props =>
    `The collection of drug sets was then visualized with a Supervenn diagram${''/* [FIGURE]*/}.`
  )
  .build()
