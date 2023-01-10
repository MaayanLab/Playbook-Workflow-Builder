
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
      <div className="flex-grow flex flex-col justify-center">
        <h2 className="text-2xl font-bold mb-2">Contributors</h2>
        <br />
        <div className="max-w-4xl flex flex-row flex-wrap self-center justify-evenly align-middle gap-4">
          <a href="https://exrna.org/" target="__blank">
            <img className="h-16" src="/logos/exRNA.png" />
          </a>
          <a href="https://www.glygen.org/" target="__blank">
            <img className="h-16" src="/logos/glygen.svg" />
          </a>
          <a href="https://lincs-dcic.org/" target="__blank">
            <img className="h-16" src="/logos/LINCS-DCIC.png" />
          </a>
          <a href="https://kidsfirstdrc.org/" target="__blank">
            <img className="h-16" src="/logos/KidsFirst.png" />
          </a>
          <a href="https://www.metabolomicsworkbench.org/" target="__blank">
            <img className="h-16" src="/logos/Metabolomics.png" />
          </a>
        </div>
      </div>
    </>
  )
}