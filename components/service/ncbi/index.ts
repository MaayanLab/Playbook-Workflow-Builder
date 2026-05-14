import { MetaNode } from '@/spec/metanode'
import { ScoredGenes } from '@/components/core/scored'
import { archs4_icon, variable_icon } from '@/icons'
import { GeneSet } from '@/components/core/set'
import { z } from 'zod'
import * as dict from '@/utils/dict'
import * as math from '@/utils/math'
import { PMCAccessionSet } from '@/components/core/set'
import { GEOAccessionTerm } from '@/components/core/term'
import { DiseaseTerm, DrugTerm, GeneTerm, PathwayTerm, PhenotypeTerm, TissueTerm } from '@/components/core/term'

const esummary_url = 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi'
const pmcconverter_url = 'https://pmc.ncbi.nlm.nih.gov/tools/idconv/api/v1/articles'


export const GEOStudyLinkedPublicationFetch = MetaNode(`GEOStudyLinkedPublicationFetch`)
  .meta({
    label: 'Fetch GEO Study Linked Publication',
    description: 'Retrieve the PMC ID(s) of publications linked to a GEO study, if any exist.',
    icon: [variable_icon],
    external: true,
  })
  .inputs({ GEOAccessionTerm })
  .output(PMCAccessionSet)
  .resolve(async (props) => {
    const gse_id = `2${'0'.repeat(11-props.inputs.GEOAccessionTerm.length)}${props.inputs.GEOAccessionTerm.replace('GSE','')}`
    const geo_req = await fetch(
      `${esummary_url}?email=&danieljbclarkemssm@gmail.comtool=playbookworkflowbuilder&id=${gse_id}&db=gds&retmode=json&idtype=acc`, {
      headers: {
        'Content-Type': 'application/json',
        'Accept': 'application/json',
      },
      method: 'GET',
      })
    if (!geo_req.ok) throw new Error('Failed to submit study accession to GEO')
    const { result : { [gse_id] : { pubmedids: ids } } } = await geo_req.json()
    const pubmed_ids = ids as string[]
    if (!pubmed_ids) throw new Error('No linked publications found')
    const pmc_req = await fetch(
      `${pmcconverter_url}?email=&danieljbclarkemssm@gmail.comtool=playbookworkflowbuilder&ids=${pubmed_ids.join(",")}&format=json`, {
      headers: {
        'Content-Type': 'application/json',
        'Accept': 'application/json',
      },
      method: 'GET',
      })
      if (!pmc_req.ok) throw new Error('Failed to retrieve PMC accessions')
    const { records } = await pmc_req.json()
    const pub_records = records as Record<string, string | number>[]
    if (!pub_records) throw new Error('No linked publications found')
    return { set: (pub_records.map((record) => record.pmcid) as string[]) }
  })
  .story(props => ({
    abstract: `The GEO study accession was used to fetch ${props.output && props.output.set.length > 1 ? 'linked publication accessions' : 'the linked publication accession'} from PMC.`,
    methods: `The GEO study accession was used to fetch ${props.output && props.output.set.length > 1 ? 'linked publication accessions' : 'the linked publication accession'} from PMC.`,
    legend: `PMC ${props.output && props.output.set.length > 1 ? 'accessions' : 'accession'} linked to ${props.inputs?.GEOAccessionTerm}` || ``,
  }))
  .build()