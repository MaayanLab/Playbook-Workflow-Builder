/**
 * End-to-end check of the MCP server (app/mcp.ts) against a running instance.
 *
 *  tsx cli/test-mcp.ts http://localhost:3000/mcp
 *
 * This creates workflows (writes to that instance's database) and resolves steps
 *  (which call external APIs): run it against a scratch server, e.g. with
 *  DATABASE_URL=sqlite:/tmp/pwb.sqlite S3_BUCKET= STORAGE_PATH=/tmp/pwb-storage
 */
import { Client } from '@modelcontextprotocol/sdk/client/index.js'
import { StreamableHTTPClientTransport } from '@modelcontextprotocol/sdk/client/streamableHttp.js'
import assert from 'node:assert/strict'

const url = process.argv[2] ?? 'http://localhost:3000/mcp'

async function connect() {
  const client = new Client({ name: 'pwb-mcp-test', version: '1.0.0' })
  await client.connect(new StreamableHTTPClientTransport(new URL(url)))
  return client
}

type Result = { result?: any, dry_run?: any, error?: { code: string, message: string, issues?: any[] }, [k: string]: any }
async function call(client: Client, name: string, args: Record<string, unknown> = {}): Promise<Result & { isError: boolean }> {
  const res = await client.callTool({ name, arguments: args })
  const structured = (res.structuredContent ?? JSON.parse((res.content as any)[0].text)) as Result
  return { ...structured, isError: Boolean(res.isError) }
}

const tests: [string, (client: Client) => Promise<void>][] = []
const test = (name: string, fn: (client: Client) => Promise<void>) => tests.push([name, fn])

test('lists tools', async (client) => {
  const { tools } = await client.listTools()
  const names = tools.map(t => t.name).sort()
  for (const name of ['describe', 'expand', 'find_steps', 'help', 'options', 'plan', 'search_published', 'set_value', 'view']) {
    assert.ok(names.includes(name), `missing tool ${name}`)
  }
})

test('parallel calls do not collide', async () => {
  const clients = await Promise.all(Array.from({ length: 8 }, connect))
  const results = await Promise.all(clients.map(c => call(c, 'options', { q: 'gene' })))
  results.forEach(r => assert.ok(!r.isError, JSON.stringify(r.error)))
  await Promise.all(clients.map(c => c.close()))
})

test('gene -> tissue expression', async (client) => {
  const start = await call(client, 'options', { q: 'gene input' })
  const input = start.result.find((o: any) => o.type === 'Input[Gene]')
  assert.ok(input?.requires_value, 'Input[Gene] requires a value')
  const gene = await call(client, 'expand', { step: { type: 'Input[Gene]', value: 'ACE2' }, include_output: true })
  assert.ok(!gene.isError, JSON.stringify(gene.error))
  assert.equal(gene.result.type, 'Term[Gene]')
  assert.equal(gene.result.status, 'done')
  const next = gene.result.next_options.find((o: any) => o.type === 'GTExTissueExpressionFromGene')
  assert.ok(next, 'next_options includes GTExTissueExpressionFromGene')
  assert.deepEqual(next.suggested_arguments, { gene: { step_id: gene.result.step_id } })
  const expr = await call(client, 'expand', { workflow_id: gene.result.workflow_id, step: { type: next.type, arguments: next.suggested_arguments } })
  assert.ok(!expr.isError, JSON.stringify(expr.error))
  assert.notEqual(expr.result.workflow_id, gene.result.workflow_id)
  assert.equal(expr.result.parent_workflow_id, gene.result.workflow_id)
  const view = await call(client, 'view', { workflow_id: expr.result.workflow_id })
  assert.equal(view.result.steps.length, 2)
  assert.equal(view.result.steps[0].value, 'ACE2')
  assert.deepEqual(view.result.steps[1].arguments, { gene: { step_id: gene.result.step_id } })
})

test('gene set input accepts several forms', async (client) => {
  for (const value of [['ACE2', 'TP53'], 'ACE2\nTP53', 'ACE2, TP53', { description: 'mine', set: ['ACE2', 'TP53'] }, JSON.stringify({ set: ['ACE2', 'TP53'] })]) {
    const res = await call(client, 'expand', { step: { type: 'Input[Set[Gene]]', value }, include_output: true })
    assert.ok(!res.isError, `${JSON.stringify(value)}: ${JSON.stringify(res.error)}`)
    assert.equal(res.result.status, 'done')
    assert.deepEqual(res.result.output.set, ['ACE2', 'TP53'])
  }
})

test('structured errors with suggestions', async (client) => {
  const missing = await call(client, 'expand', { step: { type: 'Input[Gene]' } })
  assert.ok(missing.isError)
  assert.equal(missing.error?.issues?.[0].code, 'value_required')
  assert.ok(missing.error?.issues?.[0].suggestions?.length)

  const dataType = await call(client, 'expand', { step: { type: 'Term[Gene]', value: 'ACE2' } })
  assert.equal(dataType.error?.issues?.[0].code, 'data_type_not_step_type')
  assert.ok(dataType.error?.issues?.[0].suggestions.includes('Input[Gene]'))

  const gene = await call(client, 'expand', { step: { type: 'Input[Gene]', value: 'ACE2' } })
  const bad = await call(client, 'expand', { workflow_id: gene.result.workflow_id, step: { type: 'GTExTissueExpressionFromGene', arguments: { gene: { step_id: 'nope' } } } })
  assert.equal(bad.error?.issues?.[0].code, 'invalid_step_id')
  assert.deepEqual(bad.error?.issues?.[0].suggestions, [gene.result.step_id])

  const unknownWorkflow = await call(client, 'options', { workflow_id: 'nope' })
  assert.equal(unknownWorkflow.error?.code, 'workflow_not_found')
})

test('dry run adds nothing', async (client) => {
  const res = await call(client, 'expand', { step: { type: 'Input[Gene]', value: 'ACE2' }, dry_run: true })
  assert.ok(!res.isError, JSON.stringify(res.error))
  assert.equal(res.dry_run.valid, true)
  assert.equal(res.result, undefined)
})

test('plan: multiple inputs of the same type', async (client) => {
  // GenesetsToGMT: a prompt taking multiple gene sets (& a value naming them)
  const find = await call(client, 'find_steps', { input_type: 'Set[Gene]', limit: 100 })
  const union = find.result.find((s: any) => s.type === 'GenesetsToGMT')
  assert.ok(union, `GenesetsToGMT in ${JSON.stringify(find.result.map((s: any) => s.type))}`)
  const [arg] = Object.entries(union.arguments).find(([_, a]: any) => a.multiple) as [string, any]
  const res = await call(client, 'plan', {
    steps: [
      { ref: 'a', type: 'Input[Set[Gene]]', value: ['ACE2', 'TP53'] },
      { ref: 'b', type: 'Input[Set[Gene]]', value: ['TP53', 'EGFR'] },
      { ref: 'u', type: union.type, arguments: { [arg]: { refs: ['a', 'b'] } }, value: { terms: { 0: 'set a', 1: 'set b' }, descriptions: { 0: 'first', 1: 'second' } } },
    ],
    include_output: true,
    wait_seconds: 20,
  })
  assert.ok(!res.isError, JSON.stringify(res.error))
  assert.equal(res.result.steps.length, 3)
  const last = res.result.steps[2]
  assert.equal(last.status, 'done', JSON.stringify(last))
  const view = await call(client, 'view', { workflow_id: res.result.workflow_id })
  assert.deepEqual(view.result.steps[2].arguments[arg], { step_ids: [res.result.steps[0].step_id, res.result.steps[1].step_id] })
  // options now knows two gene sets are available
  const opts = await call(client, 'options', { workflow_id: res.result.workflow_id, step_id: res.result.steps[1].step_id, q: union.type })
  const opt = opts.result.find((o: any) => o.type === union.type)
  assert.equal(opt.disabled_because, undefined)
})

test('plan is atomic', async (client) => {
  const res = await call(client, 'plan', {
    steps: [
      { ref: 'g', type: 'Input[Gene]', value: 'ACE2' },
      { type: 'GTExTissueExpressionFromGene', arguments: { gene: { ref: 'missing' } } },
    ],
  })
  assert.ok(res.isError)
  assert.equal(res.error?.issues?.[0].path, 'steps[1].arguments.gene')
})

test('set_value re-runs downstream steps', async (client) => {
  const plan = await call(client, 'plan', {
    steps: [
      { ref: 'g', type: 'Input[Gene]', value: 'ACE2' },
      { type: 'GTExTissueExpressionFromGene', arguments: { gene: { ref: 'g' } } },
    ],
  })
  assert.ok(!plan.isError, JSON.stringify(plan.error))
  const [g, e] = plan.result.steps
  const res = await call(client, 'set_value', { workflow_id: plan.result.workflow_id, step_id: g.step_id, value: 'TP53' })
  assert.ok(!res.isError, JSON.stringify(res.error))
  assert.ok(res.result.step_id_changes[g.step_id])
  assert.ok(res.result.step_id_changes[e.step_id])
  const view = await call(client, 'view', { workflow_id: res.result.workflow_id })
  assert.equal(view.result.steps[0].value, 'TP53')
})

test('describe & find_steps', async (client) => {
  const d = await call(client, 'describe', { type: 'Input[Set[Gene]]' })
  assert.equal(d.result.kind, 'step')
  assert.ok(d.result.value.json_schema)
  const t = await call(client, 'describe', { type: 'Set[Gene]' })
  assert.equal(t.result.kind, 'data')
  assert.ok(t.result.produced_by.includes('Input[Set[Gene]]'))
  const f = await call(client, 'find_steps', { q: 'enrichment', input_type: 'Set[Gene]' })
  assert.ok(f.result.length > 0)
  const bad = await call(client, 'find_steps', { output_type: 'Input[Gene]' })
  assert.deepEqual(bad.error?.issues?.[0].suggestions, ['Term[Gene]'])
})

test('search_published returns viewable workflows', async (client) => {
  const res = await call(client, 'search_published', { search: 'gene expression' })
  assert.ok(!res.isError, JSON.stringify(res.error))
  for (const p of res.result.slice(0, 3)) {
    const view = await call(client, 'view', { workflow_id: p.workflow_id })
    assert.ok(!view.isError, `${p.workflow_id}: ${JSON.stringify(view.error)}`)
  }
})

async function main() {
  const client = await connect()
  let failed = 0
  for (const [name, fn] of tests) {
    const start = Date.now()
    try {
      await fn(client)
      console.log(`ok   ${name} (${Date.now() - start}ms)`)
    } catch (e) {
      failed += 1
      console.log(`FAIL ${name}\n  ${(e as Error).message.split('\n').join('\n  ')}`)
    }
  }
  await client.close()
  console.log(failed ? `${failed} of ${tests.length} failed` : `all ${tests.length} passed`)
  process.exit(failed ? 1 : 0)
}
main()
