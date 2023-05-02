const Viz = () => null

export default function Home() {
  return (
    <div className="bp4-running-text">
      <h1 className="bp4-heading">Welcome to Signature Commons Graph</h1>
      <p className="bp4-text-large bp4-text-muted">An interactive knowledge enrichment tool for user-driven meta-graph traversal</p>
      <hr />
      <p>
        You can start <i>extending</i> your own instance graph by clicking the <svg className="inline-block" width={32} height={32} viewBox="0 0 1 1"><rect width={1} height={1} fill="lightgrey"></rect><text x={0.5} y={0.5} fontSize="0.5px" dominantBaseline="middle" textAnchor="middle" fill="black">+</text></svg> above.
      </p>
      <p>
        To completely start over rather than expanding from Home, click the Signature Commons Graph logo.
      </p>
      <figure>
        <Viz />
        <figcaption className="text-center">A global view of the meta graph</figcaption>
      </figure>
    </div>
  )
}