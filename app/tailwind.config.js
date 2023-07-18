const colors = require('tailwindcss/colors')

/** @type {import('tailwindcss').Config} */
module.exports = {
  darkMode: 'class',
  content: [
    "./pages/**/*.{js,ts,jsx,tsx}",
    "./components/**/*.{js,ts,jsx,tsx}",
    "./fragments/**/*.{js,ts,jsx,tsx}",
    "../components/**/*.{js,ts,jsx,tsx}",
  ],
  theme: {
    fontFamily: {
      sans: ['Arial', 'sans-serif'],
    },
    extend: {},
  },
  plugins: [
    require("@tailwindcss/typography"),
    require('daisyui'),
  ],
  daisyui: {
    themes: [
      {
        'light': {
          ...require("daisyui/src/colors/themes")["[data-theme=corporate]"],
          'primary': '#B3CFFF',
          'primary-content': '#2B273A',
          'secondary': '#DDDDDD',
          'primary-content': '#2B273A',
          '--border-btn': '0px',
          '--btn-text-case': '',
          '--animation-btn': '0',
          '--animation-input': '0',
        },
      },
      {
        'dark': {
          ...require("daisyui/src/colors/themes")["[data-theme=business]"],
          '--border-btn': '0px',
          '--btn-text-case': '',
          '--animation-btn': '0',
          '--animation-input': '0',
        },
      },
    ],
  },
}
