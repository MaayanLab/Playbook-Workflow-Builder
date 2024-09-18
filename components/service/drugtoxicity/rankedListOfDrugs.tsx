import { MetaNode } from '@/spec/metanode'
import { z } from 'zod'
import python from '@/utils/python'
import { Table, Cell, Column } from '@/app/components/Table'
import { drug_icon, weighted_icon } from '@/icons'
import { downloadBlob } from '@/utils/download'
import { ScoredDrugs } from '@/components/core/scored'
import { DrugCytotoxictyChembl } from './drugyCytotoxicityChembl'
import { DrugBloodBrainBarrier } from './drugBloodBrainBarrierPermeability'


export const RankedDrugToxicity = MetaNode(`[RankedDrugToxicityTable]`)
    .meta({
        label: 'Ranked Drug Toxicity',
        description: `Drugs Ranked by Cytotoxicity, Blood Brain Barrier and DrugShot`,
        icon: [drug_icon, weighted_icon],
    })
    .codec(z.array(z.object({
        drug_name: z.string(),
        confidence_zscore: z.string(),
        logBB: z.number(),
        bbb_permeable: z.string(),
        cytotoxicity_activity_id: z.number(),
        cytotoxicity_assay_description: z.string().nullable(),
        cytotoxicity_assay_type: z.string(),
        cytotoxicity_standard_type: z.string().nullable(),
        cytotoxicity_standard_units: z.string().nullable(),
        cytotoxicity_standard_value: z.number().nullable(),
        cytotoxicity_target_pref_name: z.string().nullable(),
        cytotoxicity_pchembl_value: z.number(),
    })))
    .view(rankedListTable => {
        return (
            <div className="flex flex-col gap-2">
                {/* <div className="flex flex-row gap-2">
                    <Button>Sort Drug Name</Button>
                    <Button>Sort Confidence Score</Button>
                    <Button>Sort Cytotoxicity</Button>
                </div> */}
                <Table
                    height={500}
                    cellRendererDependencies={[rankedListTable]}
                    numRows={rankedListTable.length}
                    enableGhostCells
                    enableFocusedCell
                    downloads={{
                        JSON: () => downloadBlob(new Blob([JSON.stringify(rankedListTable)], { type: 'application/json;charset=utf-8' }), 'data.json'),
                        CSV: () => downloadBlob(new Blob([
                            [
                                `drug_name,confidence_zscore,cytotoxicity_activity_id,cytotoxicity_assay_description,cytotoxicity_assay_type,cytotoxicity_standard_type,cytotoxicity_standard_units,cytotoxicity_standard_value,cytotoxicity_target_pref_name,cytotoxicity_pchembl_value,logBB,bbb_permeable`,
                                ...(rankedListTable.map((record) => [
                                    record.drug_name,
                                    record.confidence_zscore,
                                    record.cytotoxicity_activity_id,
                                    record.cytotoxicity_assay_description,
                                    record.cytotoxicity_assay_type,
                                    record.cytotoxicity_standard_type,
                                    record.cytotoxicity_standard_units,
                                    record.cytotoxicity_standard_value,
                                    record.cytotoxicity_target_pref_name,
                                    record.cytotoxicity_pchembl_value,
                                    record.logBB,
                                    record.bbb_permeable,
                                ].join(',')))
                            ].join('\n')
                        ], { type: 'text/csv;charset=utf-8' }), 'data.csv')
                    }}
                >
                    <Column
                        name="Drug Name"
                        cellRenderer={row => <Cell key={row + ''}>{rankedListTable[row].drug_name}</Cell>}
                    />
                    <Column
                        name="Confidence ZScore"
                        cellRenderer={row => <Cell key={row + ''}>{rankedListTable[row].confidence_zscore}</Cell>}
                    />
                    <Column
                        name="Cytotoxicity Activity ID"
                        cellRenderer={row => <Cell key={row + ''}>{rankedListTable[row].cytotoxicity_activity_id}</Cell>}
                    />
                    <Column
                        name="Cytotoxicity Assay Description"
                        cellRenderer={row => <Cell key={row + ''}>{rankedListTable[row].cytotoxicity_assay_description}</Cell>}
                    />
                    <Column
                        name="Cytotoxicity Assay Type"
                        cellRenderer={row => <Cell key={row + ''}>{rankedListTable[row].cytotoxicity_assay_type}</Cell>}
                    />
                    <Column
                        name="Cytotoxicity PCHEMBL Value"
                        cellRenderer={row => <Cell key={row + ''}>{rankedListTable[row].cytotoxicity_pchembl_value}</Cell>}
                    />
                    <Column
                        name="Cytotoxicity Standard Type"
                        cellRenderer={row => <Cell key={row + ''}>{rankedListTable[row].cytotoxicity_standard_type}</Cell>}
                    />
                    <Column
                        name="Cytotoxicity Standard Units"
                        cellRenderer={row => <Cell key={row + ''}>{rankedListTable[row].cytotoxicity_standard_units}</Cell>}
                    />
                    <Column
                        name="Cytotoxicity Standard Value"
                        cellRenderer={row => <Cell key={row + ''}>{rankedListTable[row].cytotoxicity_standard_value}</Cell>}
                    />
                    <Column
                        name="LogBB"
                        cellRenderer={row => <Cell key={row + ''}>{rankedListTable[row].logBB}</Cell>}
                    />
                    <Column
                        name="BBB+/BBB-"
                        cellRenderer={row => <Cell key={row + ''}>{rankedListTable[row].bbb_permeable}</Cell>}
                    />
                </Table>
            </div>
        )
    }).build()

export const RankedListDrugToxicity = MetaNode(`[RankedListOfDrugsCandidates]`)
    .meta({
        label: 'Ranked List of Drug Candidates',
        description: 'Rank List of Drug Candidates via Cytotoxicity, Blood Brain Barrier, and Confidence Score'
    })
    .inputs({ ScoredDrugs, DrugCytotoxictyChembl, DrugBloodBrainBarrier })
    .output(RankedDrugToxicity)
    .resolve(async (props) => {
        return await python(
            'components.service.drugtoxicity.produce_ranked_drug_candidates',
            { kargs: [props.inputs.ScoredDrugs, props.inputs.DrugCytotoxictyChembl, props.inputs.DrugBloodBrainBarrier] },
            message => props.notify({ type: 'info', message }),
        )
    })
    .story(props => ({
    }))
    .build()