const colors = require('tailwindcss/colors')

/** @type {import('tailwindcss').Config} */
module.exports = {
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
        playbook: {
          ...require("daisyui/src/colors/themes")["[data-theme=corporate]"],
          'primary': '#B3CFFF',
          'primary-content': '#2B273A',
          'secondary': '#DDDDDD',
          'primary-content': '#2B273A',
          '--border-btn': '0px',
          '--btn-text-case': '',
        },
      },
      'corporate',
      'business'
    ],
    darkMode: 'business',
  },
}
