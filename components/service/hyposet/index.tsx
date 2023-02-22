import React from 'react'
import { MetaNode } from '@/spec/metanode'
import { z } from 'zod'
import { supervenn_icon } from '@/icons'
import dynamic from 'next/dynamic'
import * as dict from '@/utils/dict'
import { GMT } from '@/components/data/gene_matrix_transpose'

const ReactSupervenn = dynamic(() => import('@/app/components/supervenn'), { ssr: false })

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

export const SuperVennFromGMT = MetaNode('SuperVennFromGMT')
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
  .build()
