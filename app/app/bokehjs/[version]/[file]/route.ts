import fs from 'fs'
import path from 'path'

/**
 * BokehJS's own (prebuilt) bundles, served from the installed package (see
 *  @/components/viz/bokeh/bokeh): bundling BokehJS's modules with the app breaks it in production
 *  builds, as the minifier renames its classes, by which BokehJS identifies its models.
 *
 * The version is part of the url so that the bundles can be cached for good.
 */
const files = new Set(['bokeh.min.js', 'bokeh-widgets.min.js', 'bokeh-tables.min.js', 'bokeh-gl.min.js', 'bokeh-mathjax.min.js'])

let bokehjsDir: string | undefined
function findBokehjs() {
  if (bokehjsDir) return bokehjsDir
  // the package is installed in a node_modules of this directory or one of its parents (workspaces)
  for (let dir = process.cwd(); ; dir = path.dirname(dir)) {
    const candidate = path.join(dir, 'node_modules', '@bokeh', 'bokehjs')
    if (fs.existsSync(path.join(candidate, 'package.json'))) return bokehjsDir = candidate
    if (path.dirname(dir) === dir) throw new Error('@bokeh/bokehjs is not installed')
  }
}

export async function GET(_request: Request, props: { params: Promise<{ version: string, file: string }> }) {
  const { version, file } = await props.params
  if (!files.has(file)) return new Response('Not Found', { status: 404 })
  const dir = findBokehjs()
  const installed = JSON.parse(await fs.promises.readFile(path.join(dir, 'package.json'), 'utf8')).version
  // (an immutable url only serves the version it names)
  if (version !== installed) return new Response(`BokehJS ${installed} is installed (not ${version})`, { status: 404 })
  const body = await fs.promises.readFile(path.join(dir, 'build', 'js', file))
  return new Response(body, {
    headers: {
      'Content-Type': 'text/javascript; charset=utf-8',
      'Cache-Control': 'public, max-age=31536000, immutable',
    },
  })
}
