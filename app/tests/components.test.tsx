import 'isomorphic-fetch'
import krg from '@/app/krg'

// TODO: also test rendering

describe('components', () => {
  describe('data nodes', () => {
    for (const node of krg.getDataNodes()) {
      it(node.spec, () => {
        expect(typeof node.meta.example).toBeDefined()
        expect(node.codec.encode(node.meta.example)).toBeTruthy()
      })
    }
  })
  describe('process nodes', () => {
    const Q = krg.getNextProcess('')
    const D = new Set<string>()
    while (true) {
      const component = Q.shift()
      if (component === undefined) break
      if (D.has(component.spec)) continue
      D.add(component.spec)
      if (!D.has(component.output.spec)) {
        D.add(component.output.spec)
        Q.push(...krg.getNextProcess(component.output.spec))
      }
      if ('prompt' in component) {
        it.skip(component.spec, () => {})
      } else {
        let skip = false
        const inputs: Record<string, any> = {}
        for (const key in component.inputs) {
          const value = component.inputs[key]
          const example = Array.isArray(value) ? value[0].meta.example : value.meta.example
          if (example === undefined) {
            skip = true
            break
          }
          inputs[key] = Array.isArray(value) ? [example, example] : example
        }
        if (skip) {
          it.skip(component.spec, () => {})
        } else {
          it(component.spec, async () => {
            const result = await component.resolve({ inputs })
            expect(typeof component.output.codec.encode(result)).toBe('string')
          })
        }
      }
    }
  })
})
