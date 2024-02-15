import '../next.config'
import main, { Options } from '.';

const [_node, _script, config] = process.argv
const opts = Options.parse(JSON.parse(config ?? '{}'))

main(opts).catch(e => {
  console.error(e)
  process.exit(1)
})
