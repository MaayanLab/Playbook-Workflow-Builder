import React from 'react'
import { MetaNode } from '@/spec/metanode'
import python from '@/utils/python'
import { z } from 'zod'
import { filter_icon } from '@/icons'
import { AnnData } from '../anndata'
import { Cell, Column, Table } from '@/app/components/Table'
import { FileC } from '@/components/core/file'

export const FilterAnnDataByGeneBiotype = MetaNode('FilterAnnDataByGeneBiotype')
  .meta({
    label: 'Filter AnnData Matirix by Gene Biotype',
    description: 'Selected gene biotypes to include',
    icon: [filter_icon],
  })
  .codec(z.object({
    biotypes: z.array(z.string()),
    species: z.union([z.literal('human'), z.literal('mouse')]),
  }))
  .inputs({
    anndata: AnnData,
  })
  .output(AnnData)
  .prompt(props => {
    const set = ["unknown", "tRNA", "rRNA", "snRNA", "scRNA", "snoRNA", "protein-coding", "pseudo", "transposon", "miscRNA", "ncRNA", "biological-region", "other"]
    const [biotypes, setBiotypes] = React.useState<string[]>(props.data?.biotypes || ["protein-coding"])
    const [species, setSpecies] = React.useState<'human' | 'mouse'>(props.data?.species || 'human')
    const toggle = (item: string) =>
      setBiotypes(prev =>
        prev.includes(item) ? prev.filter(b => b !== item) : [...prev, item]
      )
    return (
      <>
        <div className="flex items-center gap-3 mb-4">
          <span className="font-semibold">Species:</span>
          <div className="flex justify-center place-items-center gap-2">
            <span className="prose max-w-none hover:cursor-pointer" onClick={() => setSpecies('human')}>Human</span>
            <input
              type="checkbox"
              className="toggle"
              checked={species === 'mouse'}
              onChange={() => setSpecies(species === 'mouse' ? 'human' : 'mouse')}
            />
            <span className="prose max-w-none hover:cursor-pointer" onClick={() => setSpecies('mouse')}>Mouse</span>
          </div>
        </div>
        <Table
          height={500}
          cellRendererDependencies={[biotypes]}
          rowHeaderCellRenderer={(row) =>
            <div className="text-center block" onClick={() => toggle(set[row])}>
              <input type="checkbox" readOnly checked={biotypes.includes(set[row])} />
            </div>
          }
          numRows={set.length}
          shape={[Object.keys(biotypes).length]}
          enableGhostCells
        >
          <Column
            name={"Biotype"}
            cellRenderer={row => <Cell key={row + ''}>{set[row]}</Cell>}
          />
        </Table>
        <button className="bp5-button bp5-large" onClick={() => 
          props.submit({ biotypes,species})
          
        }>
          Submit
        </button>
      </>
    )
  })
  .resolve(async (props) => 
    await python(
      'components.data.gene_biotype_filter.gene_biotype_filter',
      { kargs: [props.inputs.anndata], kwargs: {
        species: props.data.species,
        biotypes: props.data.biotypes,
      }},
      message => props.notify({ type: 'info', message }),
  ))
  .story(props => ({
    abstract: 'Genes from an anndata matrix were filtered to include selected gene biotypes.',
    introduction: 'The AnnData Python package handles annotated data matrices in memory and on disk\\ref{doi:10.1101/2021.12.16.473007}.',
    methods: 'Gene expression count matrices are stored alongside sample metadata annotations in H5 format using AnnData\\ref{doi:10.1101/2021.12.16.473007}.',
    tableLegend: `A table showing the gene biotypes from NCBI Gene.`,
  }))
  .build()
