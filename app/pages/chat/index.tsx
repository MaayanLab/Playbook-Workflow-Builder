import React from 'react'
import dynamic from 'next/dynamic'
import Head from 'next/head'
import Link from 'next/link'

const Layout = dynamic(() => import('@/app/fragments/playbook/layout'))

export default function Chat() {
  return (
    <Layout>
      <Head><title>Chat</title></Head>
      <main className="flex-grow container mx-auto p-4 flex flex-col gap-6 justify-center place-items-center">
        <div className="prose">
          Chat Mode is still in development, we're currently testing out several methods for obtaining responses. The buttons below allow you to test these different methods.
        </div>
        <div className="flex flex-row  justify-center place-items-center gap-6">
          <Link href="/chat/gpt-prompting">
          
            <div className="tooltip tooltip-bottom" data-tip="Works through prompts to OpenAI's completions endpoint.">
              <button className="btn btn-xl btn-primary">GPT Few-shot Prompting</button>
            </div>
          </Link>
          <Link href="/chat/gpt-embedding-clf">
            <div className="tooltip tooltip-bottom" data-tip="Works by prioritizing workflow steps based on your prompt's GPT Embedding.">
              <button className="btn btn-xl btn-primary">GPT Embedding Classifier</button>
            </div>
          </Link>
          <Link href="/chat/transformer">
            <div className="tooltip tooltip-bottom" data-tip="Works by translating english text into a PWB workflow.">
            <button className="btn btn-xl btn-primary">English -&gt; Workflow Transformer</button>
            </div>
          </Link>
        </div>
      </main>
    </Layout>
  )
}
