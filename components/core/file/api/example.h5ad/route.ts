import handler from '@/utils/next-rest'
import fs from 'fs'
import path from 'path'
export const GET = handler(async (req, res) => {
  fs.createReadStream(
    path.resolve(
      process.cwd(),
      '..',
      'components',
      'core',
      'file',
      'api',
      'example.h5ad',
      'example.h5ad'
    )
  ).pipe(res)
})
