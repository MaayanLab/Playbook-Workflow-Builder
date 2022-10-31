import React from 'react'
import { MetaNode } from '@/spec/metanode'
import { z } from 'zod'
import { tissue_icon, weighted_icon } from '@/icons'

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
  .view(tissues => (
    <table>
      <thead>
        <tr>
          <th>Tissue</th>
          <th>ZScore</th>
        </tr>
      </thead>
      <tbody>
        {tissues.map(tissue =>
          <tr key={tissue.tissue}>
            <td>{tissue.tissue}</td>
            <td>{tissue.zscore.toPrecision(3)}</td>
          </tr>
        )}
      </tbody>
    </table>
  ))
  .build()
