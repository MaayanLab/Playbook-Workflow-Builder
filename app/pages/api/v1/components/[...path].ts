/**
 * Example files used to be served by component routes at these URLs, and those
 *  URLs were recorded as the source of uploads registered from the examples
 *  (`upload.url`, resolved when a DRS object is read). Keep serving them so that
 *  existing workflows built from example files keep working.
 *
 * New example uploads use the `components.file.register` procedure instead.
 */
import fs from 'fs'
import path from 'path'
import { NotFoundError, UnsupportedMethodError } from '@/spec/error'
import handler from '@/utils/next-rest'

/** legacy url path (below /api/v1/components/) => asset (below components/) */
const legacyExampleAssets: Record<string, string> = {
  'data/anndata/example.h5ad': 'data/anndata/assets/example.h5ad',
  'data/gene_count_matrix/example.tsv': 'data/gene_count_matrix/assets/example.tsv',
  'data/gene_matrix_transpose/Aging_Perturbations_from_GEO_Brain_up.gmt': 'data/gene_matrix_transpose/assets/Aging_Perturbations_from_GEO_Brain_up.gmt',
  'data/gene_signature/GTEx_aging_signature_limma.tsv': 'data/gene_signature/assets/GTEx_aging_signature_limma.tsv',
  'data/metadata_matrix/meta.tsv': 'data/metadata_matrix/assets/meta.tsv',
  'gly_gen/glycosight': 'gly_gen/glycosight/assets/01CPTAC_OVprospective_G_JHUZ_20160317_QE_r01.mzid.gz',
  'service/ctd/example_geneSetCTD.csv': 'service/ctd/assets/example_geneSetCTD.csv',
}

export default handler(async (req, res) => {
  if (req.method !== 'GET' && req.method !== 'HEAD') throw new UnsupportedMethodError(req.method)
  const segments = Array.isArray(req.query.path) ? req.query.path : [req.query.path ?? '']
  const asset = legacyExampleAssets[segments.filter(Boolean).join('/')]
  if (asset === undefined) throw new NotFoundError()
  const file = path.resolve(process.env.APP_ROOT as string, 'components', asset)
  const stat = await fs.promises.stat(file).catch(() => undefined)
  if (!stat?.isFile()) throw new NotFoundError()
  res.setHeader('Content-Length', stat.size)
  if (req.method === 'HEAD') return void res.status(200).end()
  await new Promise<void>((resolve, reject) => {
    fs.createReadStream(file).on('error', reject).pipe(res).on('finish', () => resolve())
  })
})
