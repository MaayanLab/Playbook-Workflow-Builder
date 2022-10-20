import Link from 'next/link'

export default function Header() {
  return (
    <div className="bg-primary">
      <Link href="/">
        <img
          className="mx-3 py-2 w-52 cursor-pointer"
          src="/scg-logo.svg"
        />
      </Link>
    </div>
  )
}
