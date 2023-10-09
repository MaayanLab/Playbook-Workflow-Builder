import * as tf from '@tensorflow/tfjs'
import { loadGraphModel } from '@tensorflow/tfjs-converter'

import cache from '@/utils/global_cache'
export default cache('gpt-embedding-clf', () => {
  const ctx = {
    model: null as tf.GraphModel<string | tf.io.IOHandler> | null,
    dims: null as string[] | null,
    loaded: false
  }
  Promise.all([
    loadGraphModel(`${process.env.PUBLIC_URL}/gpt-embedding-clf.js/model.json`),
    import(`@/data/gpt-embedding-clf.js/dims.json`),
  ]).then(([model, { default: dims }]) => {
    Object.assign(ctx, { model, dims, loaded: true })
  }).catch(error => {
    console.error(error);
  })
  return ctx
})
