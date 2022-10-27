import React from 'react'
import krg from '@/app/krg'
import { useRouter } from 'next/router'
import { z } from 'zod'
import dynamic from 'next/dynamic'

const Button = dynamic(() => import('@blueprintjs/core').then(({ Button }) => Button))

export default function App() {
  const router = useRouter()
  return (
    <div>
      {krg.getNextProcess('').map(proc =>
        <div key={proc.spec}>
          {Object.keys(proc.inputs).length > 0 ?
            <span>{Object.values(proc.inputs).map((i) => i.spec).join(', ')} =&gt;&nbsp;</span>
            : null}
          <Button
            onClick={async () => {
              const req = await fetch(`/api/db/fpl`, {
                method: 'POST',
                body: JSON.stringify([{
                  type: proc.spec,
                  inputs: {},
                  data: {
                    type: proc.output.spec,
                    value: proc.output.codec.encode((proc.meta as any).default)
                  },
                }])
              })
              const res = z.string().parse(await req.json())
              router.push(`/report/${res}`)
            }}
          >{proc.spec}</Button>
          <span>&nbsp; =&gt; {proc.output.spec}</span>
        </div>
      )}
    </div>
  )
}
