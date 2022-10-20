export default function Footer() {
  return (
    <div className="grid grid-flow-row md:grid-cols-4 sm:grid-cols-2 grid-cols-1 gap-4 justify-items-center items-center text-center mb-2">
      <div className="flex flex-col grid-cols-1">
        <a className="text-gray-600" href="mailto:avi.maayan@mssm.edu">Contact Us</a>
        <a className="text-gray-600" href="https://github.com/MaayanLab/signature-commons-graph/blob/master/LICENSE">Usage License</a>
      </div>
      <div className="grid-cols-1">
        <a href="https://icahn.mssm.edu/research/bioinformatics" target="_blank">
          <img className="rounded h-20" src="/icahn_cb.png" />
        </a>
      </div>
      <div className="grid-cols-1">
        <a href="https://labs.icahn.mssm.edu/maayanlab/" target="_blank">
          <img className="rounded h-20" src="/maayanlab_logo.png" />
        </a>
      </div>
      <div className="flex flex-col grid-cols-1">
        <a className="bg-gray-500 rounded-lg text-white flex p-1 m-1 hover:no-underline hover:text-white" href="https://github.com/MaayanLab/signature-commons-graph" target="_blank">
          <img className="rounded-md w-5" src="/GitHub-Mark.png" />
          &nbsp;
          View source code
        </a>
        <a className="bg-gray-500 rounded-lg text-white flex p-1 m-1 hover:no-underline hover:text-white" href="https://github.com/MaayanLab/signature-commons-graph/issues/new" target="_blank">
          <img className="rounded-md w-5" src="/GitHub-Mark.png" />
          &nbsp;
          Submit an issue
        </a>
      </div>
    </div>
  )
}