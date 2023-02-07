import usePublicUrl from '@/utils/next-public-url'
import Link from 'next/link'

export default function Header() {
  const publicUrl = usePublicUrl()

  return (
    <div className="bg-primary">
      <Link href="/">
        <img
          className="mx-3 py-2 w-52 cursor-pointer"
          src={`${publicUrl}/scg-logo.svg`}
        />
      </Link>
    </div>
  )
}
