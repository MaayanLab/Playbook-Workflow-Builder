#!/usr/bin/env ts-node
/**
 * Usage:
 * pwb MetaNodeSpec --data=data.json --inputs.key=input-key.json --output=output.json
 * 
 * E.g.
 * alias pwb="ts-node cli/pwb.ts"
 * pwb GTExTissueExpressionFromGene --inputs.gene=<(pwb Input[Gene] --data=<(echo '"ACE2"') --output=/dev/stdout) --output=/dev/stdout
 */
import krg from '@/app/krg'
import { command, run, string, array, option, multioption, subcommands } from 'cmd-ts'
import * as dict from '@/utils/dict'
import { Codec } from '@/spec/codec'
import fs from 'fs'
import { ProcessMetaNode } from '@/spec/metanode'

function readFile(file: string) {
  return new Promise<string>((resolve, reject) =>
    fs.readFile(file, { encoding: 'utf-8' }, (err, data) => {
      if (err) reject(err)
      else resolve(data)
    })
  )
}
function loadFromFile(codec: Codec<unknown>, file: string) {
  return readFile(file).then((data) => JSON.parse(data))
}
function writeFile(file: string, data: string) {
  return new Promise<void>((resolve, reject) =>
    fs.writeFile(file, data, { encoding: 'utf-8' }, (err) => {
      if (err) reject(err)
      else resolve()
    })
  )
}
function saveToFile(codec: Codec<unknown>, file: string, data: unknown) {
  return writeFile(file, JSON.stringify(data))
}
function commandFromProcessNode(proc: ProcessMetaNode) {
  return command({
    name: `pwb component ${proc.spec}`,
    args: {
      ...('codec' in proc ? {
        data: option({
          type: string,
          long: 'data',
          short: 'd',
        }),
      } : {}),
      ...dict.init(dict.items(proc.inputs).map(({ key, value }) => ({
        key: `inputs.${key}`,
        value: Array.isArray(value) ? multioption({
          type: array(string),
          long: `inputs.${key}`,
        }) : option({
          type: string,
          long: `inputs.${key}`,
        }),
      }))),
      output: option({
        type: string,
        long: 'output',
        short: 'o',
      }),
    },
    handler: async (props) => {
      const inputs: Record<string, unknown> = {}
      await Promise.all(dict.items(props).map(async ({ key, value }) => {
        const m = /^inputs\.(.+)$/.exec(key)
        if (m === null) return
        if (value === undefined) throw new Error(`Missing --${key}`)
        const input_metanode = proc.inputs[m[1]]
        if (!input_metanode) throw new Error(`Unknown --${key}`)
        if (Array.isArray(input_metanode)) {
          if (!Array.isArray(value) || value.length < 2) throw new Error(`Expected multiple --${key}`)
          inputs[m[1]] = await Promise.all(value.map(async (v) => await loadFromFile(input_metanode[0].codec, v)))
        } else {
          inputs[m[1]] = await loadFromFile(input_metanode.codec, value)
        }
      }))
      let data: unknown | undefined
      if ('codec' in proc) {
        if (props.data === undefined) throw new Error(`Missing data`)
        data = await loadFromFile(proc.codec, props.data)
      }
      const output = await proc.resolve({
        data, inputs,
        notify(status) {
          console.error(status)
        }
      })
      await saveToFile(proc.output.codec, props.output, output)
    },
  })
}

const app = subcommands({
  name: 'pwb',
  cmds: dict.init(krg.getProcessNodes().map(proc => ({
    key: proc.spec,
    value: commandFromProcessNode(proc),
  })))
})
run(app, process.argv.slice(2))
