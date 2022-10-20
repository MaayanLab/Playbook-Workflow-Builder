const colors = require('tailwindcss/colors')

module.exports = {
  theme: {
    colors: {
      primary: colors.indigo,
      secondary: colors.yellow,
      neutral: colors.gray,
    }
  }
}
/** @type {import('tailwindcss').Config} */
module.exports = {
  content: [
    "./pages/**/*.{js,ts,jsx,tsx}",
    "./components/**/*.{js,ts,jsx,tsx}",
    "./fragments/**/*.{js,ts,jsx,tsx}",
  ],
  theme: {
    colors: {
      primary: '#B3CFFF',
      secondary: '#DDDDDD',
    },
    extend: {},
  },
  plugins: [],
}
