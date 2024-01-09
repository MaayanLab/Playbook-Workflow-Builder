import handler from '@/utils/next-rest'
import fs from 'fs'
import { exampleFile } from '.'

export const GET = handler(async (req, res) => {
  fs.createReadStream(exampleFile).pipe(res)
})
