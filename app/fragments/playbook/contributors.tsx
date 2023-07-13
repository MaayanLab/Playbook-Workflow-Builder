import usePublicUrl from "@/utils/next-public-url"

export default function Contributors() {
  const publicUrl = usePublicUrl()

  return (
    <div className="flex-grow flex flex-col justify-center">
      <h2 className="text-2xl font-bold mb-2 prose">Contributors</h2>
      <br />
      <div className="max-w-4xl flex flex-row flex-wrap self-center justify-evenly align-middle gap-4">
        <a href="https://exrna.org/" target="__blank">
          <img className="h-20 p-2 rounded dark:bg-white" src={`${publicUrl}/logos/exRNA.png`} />
        </a>
        <a href="https://www.glygen.org/" target="__blank">
          <img className="h-20 p-2 rounded dark:bg-white" src={`${publicUrl}/logos/glygen.svg`} />
        </a>
        <a href="https://lincs-dcic.org/" target="__blank">
          <img className="h-20 p-2 rounded dark:bg-white" src={`${publicUrl}/logos/LINCS-DCIC.png`} />
        </a>
        <a href="https://kidsfirstdrc.org/" target="__blank">
          <img className="h-20 p-2 rounded dark:bg-white" src={`${publicUrl}/logos/KidsFirst.png`} />
        </a>
        <a href="https://www.metabolomicsworkbench.org/" target="__blank">
          <img className="h-20 p-2 rounded dark:bg-white" src={`${publicUrl}/logos/Metabolomics.png`} />
        </a>
      </div>
    </div>
  )
}