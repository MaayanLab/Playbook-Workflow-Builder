import { MetaNode } from '@/spec/metanode'
import { z } from 'zod'
import python from '@/utils/python'
import { Table, Cell, Column } from '@/app/components/Table'
import { drug_icon, weighted_icon } from '@/icons'
import { downloadBlob } from '@/utils/download'
import { ScoredDrugs } from '@/components/core/scored'


export const DrugBloodBrainBarrier = MetaNode(`[DrugBloodBrainBarrierTable]`)
    .meta({
        label: 'Drug Blood Brain Barrier Permeability',
        description: `Drug Blood Brain Barrier Permeability`,
        color: '#98D7C2',
        icon: [drug_icon, weighted_icon]
    })
    .codec(z.array(z.object({
        drug_name: z.string(),
        logBB: z.number(),
        bbb_permeable: z.string(),
    })))
    .view(drugBloodBrainBarrierTable => {
        return (
            <Table
                height={500}
                cellRendererDependencies={[drugBloodBrainBarrierTable]}
                numRows={drugBloodBrainBarrierTable.length}
                enableGhostCells
                enableFocusedCell
                downloads={{
                    JSON: () => downloadBlob(new Blob([JSON.stringify(drugBloodBrainBarrierTable)], { type: 'application/json;charset=utf-8' }), 'data.json'),
                    CSV: () => downloadBlob(new Blob([
                        [
                            `Activity ID, Assay Description, Assay Type, Drug Name, Standard Type, Standard Units, Standard Value, Target pref Name, PCHEMBL Value`,
                            ...(drugBloodBrainBarrierTable.map((record) => [
                                record.drug_name,
                                record.logBB,
                                record.bbb_permeable,
                            ].join(',')))
                        ].join('\n')
                    ], { type: 'text/csv;charset=utf-8' }), 'data.csv')
                }}
            >
                <Column
                    name="Drug Name"
                    cellRenderer={row => <Cell key={row + ''}>{drugBloodBrainBarrierTable[row].drug_name}</Cell>}
                />
                <Column
                    name="LogBB"
                    cellRenderer={row => <Cell key={row + ''}>{drugBloodBrainBarrierTable[row].logBB}</Cell>}
                />
                <Column
                    name="BBB+/BBB-"
                    cellRenderer={row => <Cell key={row + ''}>{drugBloodBrainBarrierTable[row].bbb_permeable}</Cell>}
                />
            </Table>
        )
    }).build()


export const QueryDrugBloodBrainBarrier = MetaNode('QueryDrugBloodBrainBarrier')
    .meta({
        label: 'Query Drug Blood Brain Barrier Permeability',
        description: 'Use B3DB to query Blood Brain Barrier Permeability',
        pagerank: 1,
        external: true,
    })
    .inputs({ ScoredDrugs })
    .output(DrugBloodBrainBarrier)
    .resolve(async (props) => {
        return await python(
            'components.service.drugtoxicity.query_drug_bbb_from_b3db',
            { kargs: [props.inputs.ScoredDrugs] },
            message => props.notify({ type: 'info', message }),
        )
    })
    .story(props => ({
        abstract: `Query the B3DB Database for blood-brain barrier (BBB) permeability data \\ref{doi:10.1038/s41597-021-01069-5}.`
    }))
    .build()