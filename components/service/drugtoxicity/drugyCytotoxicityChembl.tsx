import { MetaNode } from '@/spec/metanode'
import { z } from 'zod'
import python from '@/utils/python'
import { Table, Cell, Column } from '@/app/components/Table'
import { drug_icon, weighted_icon } from '@/icons'
import { downloadBlob } from '@/utils/download'
import { ScoredDrugs } from '@/components/core/scored'

export const DrugCytotoxictyChembl = MetaNode(`[DrugCytotoxictyChemblTable]`)
    .meta({
        label: 'Drug Cytotoxicty (CHEMBL)',
        description: `Drug Cytotoxicty CHEMBL`,
        color: '#98D7C2',
        icon: [drug_icon, weighted_icon]
    })
    .codec(z.array(z.object({
        activity_id: z.number(),
        assay_description: z.string().nullable(),
        assay_type: z.string(),
        molecule_pref_name: z.string(),
        standard_type: z.string().nullable(),
        standard_units: z.string().nullable(),
        standard_value: z.number().nullable(),
        target_pref_name: z.string().nullable(),
        pchembl_value: z.number(),
    })))
    .view(drugCytotoxicityTable => {
        return (
            <Table
                height={500}
                cellRendererDependencies={[drugCytotoxicityTable]}
                numRows={drugCytotoxicityTable.length}
                enableGhostCells
                enableFocusedCell
                downloads={{
                    JSON: () => downloadBlob(new Blob([JSON.stringify(drugCytotoxicityTable)], { type: 'application/json;charset=utf-8' }), 'data.json'),
                    CSV: () => downloadBlob(new Blob([
                        [
                            `Activity ID, Assay Description, Assay Type, Drug Name, Standard Type, Standard Units, Standard Value, Target pref Name, PCHEMBL Value`,
                            ...(drugCytotoxicityTable.map((record) => [
                                record.activity_id,
                                record.assay_description,
                                record.assay_type,
                                record.molecule_pref_name,
                                record.standard_type,
                                record.standard_units,
                                record.standard_value,
                                record.target_pref_name,
                                record.pchembl_value,
                            ].join(',')))
                        ].join('\n')
                    ], { type: 'text/csv;charset=utf-8' }), 'data.csv')
                }}
            >
                <Column
                    name="Drug Name"
                    cellRenderer={row => <Cell key={row + ''}>{drugCytotoxicityTable[row].molecule_pref_name}</Cell>}
                />
                <Column
                    name="Activity ID"
                    cellRenderer={row => <Cell key={row + ''}>{drugCytotoxicityTable[row].activity_id}</Cell>}
                />
                <Column
                    name="Assay Description"
                    cellRenderer={row => <Cell key={row + ''}>{drugCytotoxicityTable[row].assay_description}</Cell>}
                />
                <Column
                    name="Assay Type"
                    cellRenderer={row => <Cell key={row + ''}>{drugCytotoxicityTable[row].assay_type}</Cell>}
                />
                <Column
                    name="PCHEMBL Value"
                    cellRenderer={row => <Cell key={row + ''}>{drugCytotoxicityTable[row].pchembl_value}</Cell>}
                />
                <Column
                    name="Standard Type"
                    cellRenderer={row => <Cell key={row + ''}>{drugCytotoxicityTable[row].standard_type}</Cell>}
                />
                <Column
                    name="Standard Units"
                    cellRenderer={row => <Cell key={row + ''}>{drugCytotoxicityTable[row].standard_units}</Cell>}
                />
                <Column
                    name="Standard Value"
                    cellRenderer={row => <Cell key={row + ''}>{drugCytotoxicityTable[row].standard_value}</Cell>}
                />
            </Table>
        )
    }).build()

export const QueryDrugCytotoxocityCHEMBL = MetaNode('QueryDrugCytotoxocityCHEMBL')
    .meta({
        label: 'Query CHEMBL Drug Cytotoxocity',
        description: 'Use CHEMBL API Activity Model to obtain drug cytotoxicity data',
        pagerank: 1,
        external: true,
    })
    .inputs({ ScoredDrugs })
    .output(DrugCytotoxictyChembl)
    .resolve(async (props) => {
        return await python(
            'components.service.drugtoxicity.query_drug_cytotoxicty_from_chembl',
            { kargs: [props.inputs.ScoredDrugs] },
            message => props.notify({ type: 'info', message }),
        )
    })
    .story(props => ({
        abstract: `Cytotoxicity Results where queried for ${props.inputs?.ScoredDrugs.length} Drugs via the CHEMBL API`
    }))
    .build()