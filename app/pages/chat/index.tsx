import React from 'react'
import dynamic from 'next/dynamic'
import Head from 'next/head'
import Link from 'next/link'
import { useSessionWithId } from '@/app/extensions/next-auth/hooks'

const Layout = dynamic(() => import('@/app/fragments/playbook/layout'))

export default function Chat() {
  const { data } = useSessionWithId()
  return (
    <Layout>
      <Head><title>Chat</title></Head>
      <main className="flex-grow container mx-auto p-4 flex flex-col gap-6 justify-center place-items-center">
        <div className="prose">
          <p>Chat Mode is still in development, we're currently testing out several methods for obtaining responses. The buttons below allow you to test these different methods.</p>
        </div>
        {!data?.user?.id ?
          <Link href="/api/auth/signin?redirect=/chat">
            <button
              className="btn btn-warning"
            ><b>Please log in to use this feature</b></button>
          </Link>
          : null}
        <div className="flex flex-row flex-wrap justify-center place-items-center gap-6">
          <Link href="/chat/gpt-prompting">
            <div className="tooltip tooltip-bottom" data-tip="Works through prompts to OpenAI's completions endpoint.">
              <button className="btn btn-xl btn-primary" disabled={!data?.user?.id}>
                GPT Few-shot Prompting
              </button>
            </div>
          </Link>
          <Link href="/chat/gpt-embedding-clf">
            <div className="tooltip tooltip-bottom" data-tip="Works by prioritizing workflow steps based on your prompt's GPT Embedding.">
              <button className="btn btn-xl btn-primary" disabled={!data?.user?.id}>
                GPT Embedding Classifier
              </button>
            </div>
          </Link>
          <Link href="/chat/transformer">
            <div className="tooltip tooltip-bottom" data-tip="Works by translating english text into a PWB workflow.">
            <button className="btn btn-xl btn-primary" disabled={!data?.user?.id}>
              English -&gt; Workflow Transformer
            </button>
            </div>
          </Link>
          <Link href="/chat/transformer-embed">
            <div className="tooltip tooltip-bottom" data-tip="Works by translating english text into a PWB workflow.">
            <button className="btn btn-xl btn-primary" disabled={!data?.user?.id}>
              English + GPT Embedding -&gt; Workflow Transformer
            </button>
            </div>
          </Link>
        </div>
      </main>
    </Layout>
  )
}
