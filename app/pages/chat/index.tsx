import React from 'react'
import Head from 'next/head'
import dynamic from 'next/dynamic'
import '@webscopeio/react-textarea-autocomplete/style.css'
import krg from '@/app/krg'
import * as array from '@/utils/array'
import * as dict from '@/utils/dict'

const Layout = dynamic(() => import('@/app/fragments/playbook/layout'))
const ReactTextareaAutocomplete = dynamic(() => import('@webscopeio/react-textarea-autocomplete'))

const Item = ({ entity }: any) => <div key={entity}>{entity}</div>;

export default function Chat() {
  const [value, setValue] = React.useState('')
  const processOutputLookup = React.useMemo(() => dict.init(krg.getProcessNodes().map(proc => ({ key: proc.meta.description, value: proc.output.spec }))), [])
  return (
    <>
      <Head>
        <title>Playbook Partnership: Chat</title>
      </Head>
      <Layout>
        <main className="flex-grow container mx-auto py-4 flex flex-col">
          <ReactTextareaAutocomplete
            loadingComponent={() => <div>Loading...</div>}
            style={{
              fontSize: "18px",
              lineHeight: "20px",
              padding: 5
            }}
            containerStyle={{
              marginTop: 20,
              width: 400,
              height: 100,
              margin: "20px auto"
            }}
            minChar={0}
            trigger={{
              " ": {
                dataProvider: (token: string) => {
                  const sentences = value.split('.  ').filter(sentence => sentence in processOutputLookup)
                  const currentProc = processOutputLookup[sentences[sentences.length-1]] || ''
                  let next = krg.getNextProcess(currentProc)
                  if (next.length === 0) next = krg.getNextProcess('')
                  return array.unique(
                    next
                      .map(proc => `${proc.meta.description}.  `)
                      .filter(label => label.toLowerCase().startsWith(token.toLowerCase()))
                  )
                },
                component: Item,
                output: (item, trigger) => item as string
              }
            }}
            value={value}
            onChange={evt => {setValue(evt.target.value)}}
          />
        </main>
      </Layout>
    </>
  )
}