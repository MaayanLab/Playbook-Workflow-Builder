/**
 * This is a bit hacky but the nicer way of dealing with stories.
 * The main issue is that useMetapathInputs/useMetapathOutputs relies on SWR which only works
 *  as a React Hook. So it becomes a lot simpler to have independent processes independently
 *  assemble the story into one element -- then we can just extract the full story from the
 *  parent DOM element. This is nicer than waiting for the full story to be ready and doesn't
 *  require us to re-implement SWR.
 */
import React from 'react'
import type KRG from '@/core/KRG'
import type { Metapath } from '@/app/fragments/metapath'
import { useMetapathInputs, useMetapathOutput } from '@/app/fragments/metapath'
import extractCitations from '@/utils/citations'

/**
 * Attempt to compute the story for any given step
 */
function StoryNode({ krg, head, onChange }: { krg: KRG, head: Metapath, onChange: () => void }) {
  const processNode = krg.getProcessNode(head.process.type)
  const { data: inputs, error: inputsError } = useMetapathInputs(krg, head)
  const { data: { output, outputNode }, error: outputError } = useMetapathOutput(krg, head)
  const story = React.useMemo(() => {
    if (!processNode.story) return null
    try {
      return <>
        {processNode.story({
          inputs: !inputsError ? inputs : undefined,
          output: !outputError && outputNode.spec !== 'Error' ? output : undefined,
        })}&nbsp;
      </>
    } catch (e) {
      return null
    }
  }, [inputs, processNode, output])
  React.useEffect(onChange, [story])
  return story
}

const StoryContext = React.createContext('')

/**
 * We construct the story in a hidden element and sync changes to that element
 *  to a state variable which we expose through a context variable
 */
export function StoryProvider({ children, metapath, krg }: React.PropsWithChildren<{ metapath: Metapath[], krg: KRG }>) {
  const ref = React.useRef<HTMLSpanElement>(null)
  const [story, setStory] = React.useState('')
  React.useEffect(() => {
    setStory(() => extractCitations(ref.current?.textContent || ''))
  }, [ref, metapath])
  return (
    <>
      <StoryContext.Provider value={story}>
        {children}
      </StoryContext.Provider>
      <span ref={ref} className="hidden">
        {(metapath||[]).map(head =>
          <StoryNode key={head.id} krg={krg} head={head} onChange={() => {
            setStory(() => extractCitations(ref.current?.textContent || ''))
          }} />
        )}
      </span>
    </>
  )
}

/**
 * Any descendent of StoryProvider can access the full story text
 *  with useStory
 */
export function useStory() {
  return React.useContext(StoryContext)
}
