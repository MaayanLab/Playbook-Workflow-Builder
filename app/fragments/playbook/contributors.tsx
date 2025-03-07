import usePublicUrl from "@/utils/next-public-url"

export default function Contributors() {
  const publicUrl = usePublicUrl()

  return (
    <div className="flex-grow flex flex-col justify-center">
      <div className="max-w-4xl flex flex-row flex-wrap self-center justify-evenly align-middle gap-4">
        <a href="https://cfde.cloud/" target="_blank">
          <img className="h-20 p-2 rounded dark:bg-white" src={`${publicUrl}/logos/CFDE.png`} />
        </a>
        <a href="https://exrna.org/" target="_blank">
          <img className="h-20 p-2 rounded dark:bg-white" src={`${publicUrl}/logos/exRNA.png`} />
        </a>
        <a href="https://www.glygen.org/" target="_blank">
          <img className="h-20 p-2 rounded dark:bg-white" src={`${publicUrl}/logos/glygen.svg`} />
        </a>
        <a href="https://lincs-dcic.org/" target="_blank">
          <img className="h-20 p-2 rounded dark:bg-white" src={`${publicUrl}/logos/LINCS-DCIC.png`} />
        </a>
        <a href="https://kidsfirstdrc.org/" target="_blank">
          <img className="h-20 p-2 rounded dark:bg-white" src={`${publicUrl}/logos/KidsFirst.png`} />
        </a>
        <a href="https://www.metabolomicsworkbench.org/" target="_blank">
          <img className="h-20 p-2 rounded dark:bg-white" src={`${publicUrl}/logos/Metabolomics.png`} />
        </a>
      </div>
    </div>
  )
}