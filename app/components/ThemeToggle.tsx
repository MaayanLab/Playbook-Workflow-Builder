import React from 'react'
import useLocalStorage from '@/utils/localstorage'

export default function ThemeToggle() {
  const [theme, setTheme] = useLocalStorage('light' as 'light' | 'dark')
  React.useEffect(() => {
    if (theme === 'dark') {
      document.getElementsByTagName('html')[0].setAttribute('data-theme', 'dark')
      document.getElementsByTagName('html')[0].classList.add('dark')
      document.getElementsByTagName('html')[0].classList.add('bp4-dark')
      return () => {
        document.getElementsByTagName('html')[0].setAttribute('data-theme', 'light')
        document.getElementsByTagName('html')[0].classList.remove('dark')
        document.getElementsByTagName('html')[0].classList.remove('bp4-dark')
      }
    }
  }, [theme])
  return (
    <input
      type="checkbox"
      className="toggle"
      title="Toggle light/dark theme"
      checked={theme === 'dark'}
      onClick={evt => {setTheme(evt.currentTarget.checked ? 'dark' : 'light')}}
    />
  )
}
