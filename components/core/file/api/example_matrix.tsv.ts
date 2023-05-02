import { UnsupportedMethodError } from '@/spec/error'
import handler from '@/utils/next-rest'
import fs from 'fs'
export default handler(async (req, res) => {
  if (req.method !== 'GET') throw new UnsupportedMethodError()
  fs.createReadStream('../components/core/file/public/example_matrix.tsv')
    .pipe(res)
})
