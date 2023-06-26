import handler from '@/utils/next-rest'
import fs from 'fs'
import path from 'path'
export const GET = handler(async (req, res) => {
  fs.createReadStream(
    path.resolve(
      process.env.APP_ROOT as string,
      'components',
      'core',
      'file',
      'api',
      'example_matrix.tsv',
      'example_matrix.tsv'
    )
  ).pipe(res)
})
