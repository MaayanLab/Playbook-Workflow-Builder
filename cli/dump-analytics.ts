import db from '@/app/db'
import fpprg from '@/app/fpprg'

async function main() {
  (await Promise.all(
    (await db.objects.fpl.findMany())
      .map(async (fpl_record) => {
        const fpl_object = await fpprg.getFPL(fpl_record.id)
        if (fpl_object) return { fpl_record, fpl_object }
      })
  )).forEach((fpl) => {
    if (fpl) {
      const { fpl_record, fpl_object } = fpl
      console.log([
        fpl_record.hits,
        ...fpl_object.resolve().map(fpl => fpl.process.type)
      ].join('\t'))
    }
  })
}
main()
  .then(() => {process.exit(0)})
  .catch((err) => {console.error(err)})
