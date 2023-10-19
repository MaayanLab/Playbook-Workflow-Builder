import React from 'react'
import useLocalStorage from '@/utils/localstorage'

export default function ThemeToggle() {
  const [theme, setTheme] = useLocalStorage('light' as 'light' | 'dark')
  React.useEffect(() => {
    if (theme === 'dark') {
      document.getElementsByTagName('html')[0].setAttribute('data-theme', 'dark')
      document.getElementsByTagName('html')[0].classList.add('dark')
      document.getElementsByTagName('html')[0].classList.add('bp5-dark')
      return () => {
        document.getElementsByTagName('html')[0].setAttribute('data-theme', 'light')
        document.getElementsByTagName('html')[0].classList.remove('dark')
        document.getElementsByTagName('html')[0].classList.remove('bp5-dark')
      }
    }
  }, [theme])
  return (
    <input
      type="checkbox"
      className="toggle"
      title="Toggle light/dark theme"
      checked={theme === 'dark'}
      onChange={evt => {setTheme(evt.currentTarget.checked ? 'dark' : 'light')}}
    />
  )
}
