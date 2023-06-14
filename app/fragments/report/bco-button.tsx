import React from 'react'
import dynamic from 'next/dynamic'
import useSWRMutation from 'swr/mutation'
import { biocompute_icon } from '@/icons'

const Bp4Popover = dynamic(() => import('@blueprintjs/popover2').then(({ Popover2 }) => Popover2))
const Bp4Menu = dynamic(() => import('@blueprintjs/core').then(({ Menu }) => Menu))
const Bp4MenuItem = dynamic(() => import('@blueprintjs/core').then(({ MenuItem }) => MenuItem))
const Bp4Spinner = dynamic(() => import('@blueprintjs/core').then(({ Spinner }) => Spinner))
const Bp4Alert = dynamic(() => import('@blueprintjs/core').then(({ Alert }) => Alert))
const Icon = dynamic(() => import('@/app/components/icon'))

async function submitBcoDraft(endpoint: string): Promise<{ object_id: string }> {
  const req = await fetch(endpoint, { method: 'POST' })
  if (req.status !== 200) throw new Error(await req.json())
  else return await req.json()
}

export default function BCOButton({ session_id, id, metadata, disabled }: { session_id?: string, id?: string, metadata: { title: string, description: string | undefined }, disabled: boolean }) {
  const { trigger, isMutating, error } = useSWRMutation(id ? `${session_id ? `/api/socket/${session_id}` : ''}/api/bco/${id}/draft` : null, submitBcoDraft)
  const [isError, setIsError] = React.useState(false)
  React.useEffect(() => {
    if (error) {
      setIsError(() => true)
      console.error(error)
    }
  }, [error])
  if (!id) return null
  else if (isMutating) return <Bp4Spinner className="inline-block" size={20} />
  return (
    <>
      <Bp4Popover
        className={disabled ? 'cursor-not-allowed' : 'cursor-pointer'}
        disabled={disabled}
        content={
          <Bp4Menu>
            <a href={`/api/bco/${id}?metadata=${encodeURIComponent(JSON.stringify(metadata))}`} download={`${metadata.title.replace(/ /g, '-')}-bco.json`}>
              <Bp4MenuItem
                icon="document"
                text="Download BCO"
              />
            </a>
            <Bp4MenuItem
              icon="send-to"
              text="Draft in BioCompute Portal"
              onClick={async (evt) => {
                trigger().then((res) => {
                  if (res) window.open(`https://biocomputeobject.org/builder?${res.object_id}`, '_blank')
                })
              }}
            />
          </Bp4Menu>
        }
        placement="bottom"
      >
        <Icon
          icon={biocompute_icon}
          className={disabled ? 'fill-gray-400' : 'fill-black dark:fill-white'}
          title={disabled ? 'Save to Create BCO' : 'Create BCO'}
        />
      </Bp4Popover>
      <Bp4Alert
        confirmButtonText="Okay"
        icon="error"
        intent="danger"
        isOpen={isError}
        canEscapeKeyCancel
        canOutsideClickCancel
        onCancel={() => {setIsError(false)}}
        onConfirm={() => {setIsError(false)}}
      >
        An error ocurred while attempting to register the BCO.
        Please try again later.
      </Bp4Alert>
    </>
  )
}
