/**
 * Link a step which can't run to the steps of the report which keep it from running: the inputs
 *  in error (for "... can't run because of error in ..."), or the prompts waiting for input.
 */
import React from 'react'
import type KRG from '@/core/KRG'
import { type Metapath, useFPL, useMetapathOutput } from '@/app/fragments/metapath'
import * as dict from '@/utils/dict'
import { ExLink } from '@/app/fragments/ex-router'

/**
 * Where a link to a step goes: in a report, its cell (on the same page); in the graph view, which
 *  shows one step at a time, its node
 */
export type StepLinkMode = { mode: 'report' } | { mode: 'graph', session_id?: string }

/** The error engine.ts gives a step when one of its inputs is in error */
const INPUT_ERROR = /^(?:Error: )?([\s\S]*?) can't run because of error in ([\s\S]+)$/

export function isInputError(output: unknown): output is string {
  return typeof output === 'string' && INPUT_ERROR.test(output)
}

function useSteps(id: string) {
  const { data: metapath } = useFPL(id)
  return React.useMemo(() => new Map((metapath ?? []).map(step => [step.process.id, step])), [metapath])
}

function inputsOf(steps: Map<string, Metapath>, step: Metapath) {
  return dict.sortedValues(step.process.inputs).flatMap(({ id }) => steps.has(id) ? [steps.get(id) as Metapath] : [])
}

function StepLinks({ krg, id, steps, link }: { krg: KRG, id: string, steps: Metapath[], link: StepLinkMode }) {
  return <>{steps.map((step, i) => {
    const label = krg.getProcessNode(step.process.type)?.meta.label ?? step.process.type
    return (
      <React.Fragment key={step.id}>
        {i > 0 ? ', ' : null}
        {link.mode === 'graph'
          ? <ExLink className="link" href={`${link.session_id ? `/session/${link.session_id}` : ''}/graph/${id}/node/${step.id}`} shallow>{label}</ExLink>
          : <a className="link" href={`#${step.id}:data`}>{label}</a>}
      </React.Fragment>
    )
  })}</>
}

/** Reports whether a step's output is an error (rendering nothing) */
function ErrorProbe({ krg, step, onChange }: { krg: KRG, step: Metapath, onChange: (id: string, error: boolean) => void }) {
  const { data: { outputNode } } = useMetapathOutput({ krg, head: step })
  const error = outputNode?.spec === 'Error'
  React.useEffect(() => { onChange(step.id, error) }, [step.id, error, onChange])
  return null
}

/** "X can't run because of error in Y", with Y linked to the input step(s) in error */
export function InputError({ krg, id, head, output, link = { mode: 'report' } }: { krg: KRG, id: string, head: Metapath, output: string, link?: StepLinkMode }) {
  const [, what, type] = INPUT_ERROR.exec(output) ?? [undefined, output, '']
  const steps = useSteps(id)
  const inputs = React.useMemo(() => inputsOf(steps, head), [steps, head])
  const [errors, setErrors] = React.useState<Record<string, boolean>>({})
  const onChange = React.useCallback((id: string, error: boolean) => setErrors(errors => errors[id] === error ? errors : { ...errors, [id]: error }), [])
  const failed = inputs.filter(step => errors[step.id])
  return (
    <div className="alert alert-error">
      {inputs.map(step => <ErrorProbe key={step.id} krg={krg} step={step} onChange={onChange} />)}
      <span>{what} can't run because of an error in {failed.length ? <StepLinks krg={krg} id={id} steps={failed} link={link} /> : type}</span>
    </div>
  )
}

/** The nearest prompts upstream of a step which don't have a value yet */
function awaitingInput(steps: Map<string, Metapath>, krg: KRG, head: Metapath) {
  const found = new Map<string, Metapath>()
  const seen = new Set<string>()
  const visit = (step: Metapath) => {
    for (const input of inputsOf(steps, step)) {
      if (seen.has(input.id)) continue
      seen.add(input.id)
      const node = krg.getProcessNode(input.process.type)
      if (node && 'prompt' in node && !input.process.data) found.set(input.id, input)
      else visit(input)
    }
  }
  visit(head)
  return [...found.values()]
}

/** "Waiting for input", from the prompt(s) which need it */
export function WaitingForInput({ krg, id, head, link = { mode: 'report' } }: { krg: KRG, id: string, head: Metapath, link?: StepLinkMode }) {
  const steps = useSteps(id)
  const waiting = React.useMemo(() => awaitingInput(steps, krg, head), [steps, krg, head])
  if (!waiting.length) return <>Waiting for input</>
  return <>Waiting for input in <StepLinks krg={krg} id={id} steps={waiting} link={link} /></>
}
