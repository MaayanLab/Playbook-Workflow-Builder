export default function Footer() {
  return (
    <div className="grid grid-flow-row md:grid-cols-3 sm:grid-cols-2 grid-cols-1 gap-4 justify-items-center items-center text-center mb-2">
      <div className="flex flex-col grid-cols-1">
        <a className="text-gray-600" href="mailto:avi.maayan@mssm.edu">Contact Us</a>
        <a className="text-gray-600" href="https://github.com/nih-cfde/playbook-partnership/blob/master/LICENSE">Usage License</a>
      </div>
      <div className="grid-cols-1">
        <a href="https://www.nih-cfde.org/" target="_blank">
          <img className="rounded h-20" src="/logos/CFDE.png" />
        </a>
      </div>
      <div className="flex flex-col grid-cols-1 gap-1">
        <a className="btn btn-xs btn-secondary rounded-lg gap-1" href="https://github.com/nih-cfde/playbook-partnership" target="_blank">
          <img className="rounded-md w-4 justify-self-start" src="/GitHub-Mark.png" />
          <span className="flex-grow">View source code</span>
        </a>
        <a className="btn btn-xs btn-secondary rounded-lg gap-1" href="https://github.com/nih-cfde/playbook-partnership/issues/new" target="_blank">
          <img className="rounded-md w-4 justify-self-start" src="/GitHub-Mark.png" />
          <span className="flex-grow">Submit an issue</span>
        </a>
      </div>
    </div>
  )
}