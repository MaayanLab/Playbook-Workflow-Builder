import usePublicUrl from '@/utils/next-public-url'
import { ExLink } from '@/app/fragments/ex-router'

export default function Header() {
  const publicUrl = usePublicUrl()

  return (
    <div className="bg-primary">
      <ExLink href="/">
        <img
          className="mx-3 py-2 w-52 cursor-pointer"
          src={`${publicUrl}/scg-logo.svg`}
        />
      </ExLink>
    </div>
  )
}
