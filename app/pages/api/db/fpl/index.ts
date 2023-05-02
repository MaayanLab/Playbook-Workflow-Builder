import fpprg from '@/app/fpprg'
import { FPL, } from '@/core/FPPRG'
import { z } from 'zod'
import { IdOrProcessC } from '@/core/FPPRG'
import { NotFoundError, UnsupportedMethodError } from '@/spec/error'
import handler from '@/utils/next-rest'

const BodyType = z.array(IdOrProcessC)

export default handler(async (req, res) => {
  if (req.method !== 'POST') throw new UnsupportedMethodError()
  const processArray = await Promise.all(BodyType.parse(JSON.parse(req.body)).map(fpprg.resolveProcess))
  const processArrayFPL = FPL.fromProcessArray(processArray)
  if (!processArrayFPL) throw new NotFoundError()
  const fpl = await fpprg.upsertFPL(processArrayFPL)
  res.status(200).json(fpl.id)
})
