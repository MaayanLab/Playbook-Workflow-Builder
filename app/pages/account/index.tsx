import * as Auth from 'next-auth/react'

export default function Account() {
  const { data: session } = Auth.useSession()
  if (session && session.user) {
    return <>
      Signed in as {session.user.email} <br/>
      <button onClick={() => Auth.signOut()}>Sign out</button>
    </>
  }
  return <>
    Not signed in <br/>
    <button onClick={() => Auth.signIn()}>Sign in</button>
  </>
}
