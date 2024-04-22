import React from 'react'
import { MetaNode } from '@/spec/metanode'
import { z } from 'zod'
import python from '@/utils/python'
import { GeneTerm } from '@/components/core/term'
import { GeneInfo, GeneInfoFromGeneTerm } from '@/components/service/mygeneinfo'
import { Table, Cell, Column } from '@/app/components/Table'
import { tumor_icon } from '@/icons'
import { downloadBlob } from '@/utils/download'

export const TumorGeneExpression = MetaNode(`[TumorGeneExpression]`)
  .meta({
    label: 'Tumor Gene Expression',
    description: `Gene expression in tumors`,
    color: '#98D7C2',
    icon: [tumor_icon],
  })
  .codec(z.array(z.object({ TPM_mean : z.number(),
                            TPM_sd : z.number(),
                            TPM_median : z.number(),
                            Disease : z.string(),
                            Gene_symbol : z.string(),
                            Gene_Ensembl_ID : z.string(),
                            Dataset : z.string(),
                            zscore: z.number()})))
  .view(expressionTable => {
    return (
      <Table
        height={500}
        cellRendererDependencies={[expressionTable]}
        numRows={expressionTable.length}
        enableGhostCells
        enableFocusedCell
        downloads={{
          JSON: () => downloadBlob(new Blob([JSON.stringify(expressionTable)], { type: 'application/json;charset=utf-8' }), 'data.json'),
          CSV: () => downloadBlob(new Blob([
            [
              `Gene ID,Data Set,Disease,TPM Mean,TPM Stand. Dev.,TPM Median,Ensembl ID,Zscore`,
              ...(expressionTable.map((record) => [
                record.Gene_symbol,
                record.Dataset,
                record.Disease,
                record.TPM_mean,
                record.TPM_sd,
                record.TPM_median,
                record.Gene_Ensembl_ID,
                record.zscore,
              ].join(',')))
            ].join('\n')
          ], { type: 'text/csv;charset=utf-8' }), 'data.csv'),
        }}
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
        <Column
          name="Zscore"
          cellRenderer={row => <Cell key={row+''}>{expressionTable[row].zscore}</Cell>}
        />
      </Table>
    )
  })
  .build()

export const KFTumorExpression = MetaNode('KFTumorExpression')
  .meta({
    label: 'Query KF Gene Expression in Tumor',
    description: 'Use KF API to obtain tumors expressing the given gene',
    pagerank: 1,
  })
  .inputs({ gene_info: GeneInfo })
  .output(TumorGeneExpression)
  .resolve(async (props) => {
    return await python(
      'components.service.kf.main',
      { kargs: [props.inputs.gene_info.ensembl?.gene]},
      message => props.notify({ type: 'info', message }),
    )
  })
  .story(props => `Gene expression in tumors for ${props.inputs ? props.inputs.gene_info.symbol : 'the gene'} were queried from the Open Pediatric Cancer Atlas API [\\ref{doi:10.1016/j.xgen.2023.100340}].`)
  .build()

export const KFTumorExpressionFromGene = MetaNode('KFTumorExpressionFromGene')
  .meta(KFTumorExpression.meta)
  .inputs({ gene: GeneTerm })
  .output(KFTumorExpression.output)
  .resolve(async (props) => {
    const gene_info = await GeneInfoFromGeneTerm.resolve(props)
    return await KFTumorExpression.resolve({ ...props, inputs: { gene_info } })
  })
  .story(props => `Gene expression in tumors for ${props.inputs ? props.inputs.gene : 'the gene'} were queried from the Open Pediatric Cancer Atlas API [\\ref{doi:10.1016/j.xgen.2023.100340}].`)
  .build()
