import dynamic from 'next/dynamic'

const Contributors = dynamic(() => import('@/app/fragments/playbook/contributors'))

export default function Home() {
  return (
    <>
      <div className="bp4-running-text">
        <h2 className="bp4-heading">Construct a Playbook</h2>
        <p className="bp4-text-large bp4-text-muted">This interactive tool enables dynamic exploration and construction of user-driven workflows leveraging curated APIs.</p>
        <hr />
        <p>
          You can start <i>extending</i> your own instance graph by clicking the <svg className="inline-block" width={32} height={32} viewBox="0 0 1 1"><rect width={1} height={1} fill="lightgrey"></rect><text x={0.5} y={0.5} fontSize="0.5px" dominantBaseline="middle" textAnchor="middle" fill="black">+</text></svg> above.
        </p>
        <p>
          To completely start over rather than expanding from Home, click the logo.
        </p>
      </div>
      <Contributors />
    </>
  )
}