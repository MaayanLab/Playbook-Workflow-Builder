/**
 * Run python scripts/functions from typescript
 * Usage:
 * ```ts
 * import pythonHelper from '@/utils/python_helper'
 * async function run() {
 *   try {
 *     const result = await pythonHelper('@/path/from/root/myscript.py', 'myfunc', { kargs: [ 'myarg' ], kwargs: { 'myarg2': '2' } })
 *   } catch (e) {
 *     // e is a ProcessError
 *   }
 * }
 */
let spawn, path
if (!process.browser) {
  spawn = require('child_process').spawn
  path = require('path')
}

export class ProcessError extends Error {
  constructor(public message: string, public exitcode: number | null) {
    super(message)
    this.name = 'ProcessError'
  }
}

/**
 * @param script: of the form python_file.py:func_in_file
 * @param args: The kwargs to provide to the file
 */
export default function python<T>(script_path: string, func_name: string, args: { kargs?: unknown[], kwargs?: Record<string, unknown> }): Promise<T> {
  return new Promise((resolve, reject) => {
    let stdin: string
    try {
      stdin = JSON.stringify(args)
    } catch (e) {
      throw new ProcessError(`Process input could not be serialized`, null)
    }
    const proc = spawn('python3', [
      path.join(process.env.PYTHON_ROOT, 'utils/helper.py'),
      path.join(process.env.PYTHON_ROOT, path.relative('@', script_path)),
      func_name,
    ])
    let stdout = ''
    let stderr = ''
    proc.stdout.on('data', (chunk: string) => { stdout += chunk })
    proc.stderr.on('data', (chunk: string) => { stderr += chunk })
    proc.on('close', (code) => {
      if (code !== 0) {
        reject(new ProcessError(stderr || `Process exited with unexpected code ${code}`, code))
      } else {
        if (stderr) {
          console.warn(`${script_path}[${func_name}]: ${stderr}`)
        }
        try {
          resolve(JSON.parse(stdout))
        } catch (e) {
          reject(new ProcessError(`Process output could not be parsed as json`, code))
        }
      }
    })
    proc.stdin.end(stdin)
  })
}
