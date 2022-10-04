import React from 'react'
import { MetaNode } from '@/spec/metanode'

export type SignificantTissue = {
  /**
   * The tissue term
   */
  tissue: string,
  /**
   * A zscore fir the significance
   */
  zscore: number,
}

export const SignificantTissues = MetaNode.createData('SignificantTissues')
  .meta({
    label: 'SignificantTissues',
    description: 'Tissues scored using a combined stouffer statistic',
  })
  .codec<Array<SignificantTissue>>()
  .view(tissues => (
    <table>
      <tr>
        <th>Tissue</th>
        <th>ZScore</th>
      </tr>
      {tissues.map(tissue =>
        <tr key={tissue.tissue}>
          <td>{tissue.tissue}</td>
          <td>{tissue.zscore.toPrecision(3)}</td>
        </tr>
      )}
    </table>
  ))
  .build()
