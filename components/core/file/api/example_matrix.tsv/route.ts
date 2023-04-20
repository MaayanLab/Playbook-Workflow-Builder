import handler from '@/utils/next-rest'
import fs from 'fs'
export const GET = handler(async (req, res) => {
  fs.createReadStream('./example_matrix.tsv')
    .pipe(res)
})
