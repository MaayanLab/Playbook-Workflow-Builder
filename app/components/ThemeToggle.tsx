import React from 'react'
import useLocalStorage from '@/utils/localstorage'

export default function ThemeToggle() {
  const [theme, setTheme] = useLocalStorage('theme')
  const [currentTheme, setCurrentTheme] = React.useState<string>(() => window.matchMedia && window.matchMedia('(prefers-color-scheme: dark)').matches ? 'dark' : 'light')
  React.useEffect(() => {
    if (!window.matchMedia) return
    const media = window.matchMedia('(prefers-color-scheme: dark)')
    const listener = (evt: MediaQueryListEvent) => {
      setTheme(null)
      setCurrentTheme(evt.matches ? 'dark' : 'light')
    }
    media.addEventListener('change', listener)
    return () => {
      media.removeEventListener('change', listener)
    }
  }, [])
  React.useEffect(() => {
    if (theme === null) return
    setCurrentTheme(theme === 'dark' ? 'dark' : 'light')
  }, [theme])
  React.useEffect(() => {
    if (currentTheme === 'dark') {
      document.getElementsByTagName('html')[0].setAttribute('data-theme', 'dark')
      document.getElementsByTagName('html')[0].classList.add('dark')
      document.getElementsByTagName('html')[0].classList.add('bp5-dark')
    } else {
      document.getElementsByTagName('html')[0].setAttribute('data-theme', 'light')
      document.getElementsByTagName('html')[0].classList.remove('dark')
      document.getElementsByTagName('html')[0].classList.remove('bp5-dark')
    }
  }, [currentTheme])
  return (
    <input
      type="checkbox"
      className="toggle"
      title="Toggle light/dark theme"
      checked={currentTheme === 'dark'}
      onChange={evt => {setTheme(evt.currentTarget.checked ? 'dark' : 'light')}}
    />
  )
}
