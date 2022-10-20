import Link from 'next/link'

export default function Header() {
  return (
    <div className="bg-primary">
      <Link href="/">
        <h1 className="text-4xl font-bold p-2 cursor-pointer">Playbook Partnership Interactive Workflow Builder</h1>
      </Link>
    </div>
  )
}
