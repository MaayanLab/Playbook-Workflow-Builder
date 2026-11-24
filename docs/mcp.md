# MCP server

PWB serves a [Model Context Protocol](https://modelcontextprotocol.io) server at `/mcp`
(streamable HTTP, stateless, JSON responses) so agents can find, build, and run workflows.
The in-app chat (`app/trpc/routers/chat.ts`) uses it, and any MCP client can connect to
`https://playbook-workflow-builder.cloud/mcp`.

Source: [`app/mcp.ts`](../app/mcp.ts) (tools) and [`app/mcp-util.ts`](../app/mcp-util.ts) (helpers).
The same guide is available to agents through the `help` tool.

## Concepts

- A **workflow** is an ordered list of steps, identified by a `workflow_id`. Workflows are
  persistent: every change makes a **new** `workflow_id`, and old ids keep pointing at the
  workflow as it was. Building on an old id starts a separate branch.
- A **step** has a **step type** (a process, e.g. `Input[Gene]`, `GTExTissueExpressionFromGene`)
  and produces a value of a **data type** (e.g. `Term[Gene]`, `Scored[Tissue]`). It is identified by a `step_id`.
- Step **arguments** take the `step_id` of a prior step that produces the required data type.
  Arguments with `multiple: true` take a list: `{"step_ids": [...]}`.
- Steps with `requires_value` (input steps & other prompts) need a `value`. `describe(step_type)`
  gives its JSON schema and an example. Common forms are accepted: `"ACE2"` for a term;
  `["ACE2","TP53"]`, `"ACE2\nTP53"`, `"ACE2, TP53"` or `{"description": "...", "set": [...]}` for a set;
  JSON strings of any of these.
- Steps resolve asynchronously. `status` is `pending`, `done`, `error` or `waiting_for_input`.
  `include_output: true` waits (up to `wait_seconds`, default 10, max 60) and returns a
  summary of the output (long lists and strings are truncated with a `…(+N more)` marker).

## Tools

| tool | purpose |
| --- | --- |
| `help` | this guide |
| `search_published(search)` | published workflows. To use one: `view` it (which opens it for the user), then `set_value` its input steps; an empty result comes with a `hint` on how to search |
| `options(workflow_id?, step_id?, q?, limit?, include_disabled?)` | what can be added after the latest (or given) step; input steps when there's no workflow. Each option has `suggested_arguments` (ready for `expand`), `candidates` per argument, `requires_value`/`value_example`, or `disabled_because` |
| `find_steps(q?, input_type?, output_type?, limit?)` | search every step type, e.g. everything producing `Set[Gene]` |
| `describe(type?, q?)` | a step type (arguments, output, value format) or data type (example, schema, producers, consumers); no type lists data types |
| `expand(workflow_id?, step, dry_run?, include_output?, wait_seconds?)` | add one step; returns the new `workflow_id`, `parent_workflow_id`, status, `link` & `next_options` |
| `plan(workflow_id?, steps, dry_run?, include_output?, wait_seconds?)` | add several steps atomically; steps refer to earlier ones with `{"ref": "..."}` (refs default to the index, `"0"`, `"1"`, ...) |
| `set_value(workflow_id, step_id, value)` | change an input's value, recomputing what depends on it; returns `step_id_changes` (old → new) |
| `view(workflow_id, step_id?, include_output?, wait_seconds?)` | the steps of a workflow with arguments, values (with `requires_value` & `value_format` for `set_value`), status (and outputs). `summary.failures` lists the steps which failed themselves; other errors carry `blocked_by`, the failed step causing them |

Example: a gene's tissue expression.

```json
{"name": "expand", "arguments": {"step": {"type": "Input[Gene]", "value": "ACE2"}}}
→ {"result": {"workflow_id": "W1", "step_id": "S1", "type": "Term[Gene]", "status": "done",
   "next_options": [{"type": "GTExTissueExpressionFromGene", "suggested_arguments": {"gene": {"step_id": "S1"}}}, ...]}}

{"name": "expand", "arguments": {"workflow_id": "W1", "include_output": true,
  "step": {"type": "GTExTissueExpressionFromGene", "arguments": {"gene": {"step_id": "S1"}}}}}
→ {"result": {"workflow_id": "W2", "parent_workflow_id": "W1", "step_id": "S2", "type": "Scored[Tissue]",
   "status": "done", "output": [{"term": "Testis", "zscore": 6.06}, ..., "…(+44 more, 54 total)"], "link": ".../report/W2"}}
```

The same in one call:

```json
{"name": "plan", "arguments": {"steps": [
  {"ref": "g", "type": "Input[Gene]", "value": "ACE2"},
  {"type": "GTExTissueExpressionFromGene", "arguments": {"gene": {"ref": "g"}}}
]}}
```

## Errors

Failed calls return `isError: true` and:

```json
{"error": {
  "code": "invalid_step",
  "message": "1 problem with the requested step, nothing was added",
  "issues": [{
    "path": "step.arguments.gene", "code": "invalid_step_id",
    "message": "\"abc\" is not a step in this workflow; steps producing Term[Gene]: S1",
    "expected": "Term[Gene]", "suggestions": ["S1"]
  }],
  "hint": "Fix the issues (see suggestions) and retry. Use dry_run: true to check without adding anything."
}}
```

Issue codes include `unknown_step_type`, `data_type_not_step_type` (a data type was given as
`step.type`; suggestions are the steps producing it), `unexpected_argument`, `missing_argument`,
`invalid_step_id`, `type_mismatch`, `expected_single`, `expected_multiple`, `value_required`,
`invalid_value` (with `expected` JSON schema and an example), `unknown_data_type`. Every issue
is reported at once, and `plan` adds nothing unless all steps are valid.

## Testing

`cli/test-mcp.ts` checks a running server end to end. It creates workflows and resolves
steps (calling external APIs), so point it at a scratch server:

```bash
DATABASE_URL=sqlite:/tmp/pwb.sqlite STORAGE_PATH=/tmp/pwb-storage S3_BUCKET= S3_ENDPOINT= \
  PUBLIC_URL=http://localhost:3123 PORT=3123 npm run dev --workspace=app
npx tsx cli/test-mcp.ts http://localhost:3123/mcp
```
