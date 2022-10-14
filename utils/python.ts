/**
 * Run python scripts/functions from typescript
 * Usage:
 * ```ts
 * import python from '@/utils/python'
 * async function run() {
 *   try {
 *     const result = await python('@/path/from/root/myscript.py', 'myfunc', { kargs: [ 'myarg' ], kwargs: { 'myarg2': '2' } })
 *   } catch (e) {
 *     // e is a ProcessError
 *   }
 * }
 */
import type * as child_process_type from 'child_process'
const spawn = typeof window === 'undefined' ? (require('child_process') as typeof child_process_type).spawn : undefined
import type * as path_type from 'path'
const path = typeof window === 'undefined' ? require('path') as typeof path_type : undefined

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
export default function python<T>(pathspec: string, args: { kargs?: unknown[], kwargs?: Record<string, unknown> }): Promise<T> {
  if (typeof window !== 'undefined') throw new Error("python is server side only")
  return new Promise((resolve, reject) => {
    let stdin: string
    try {
      stdin = JSON.stringify(args)
    } catch (e) {
      throw new ProcessError(`Process input could not be serialized`, null)
    }
    const proc = spawn('python3', [
      path.join(process.env.PYTHON_ROOT, 'utils/helper.py'),
      pathspec,
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
          console.warn(`[${pathspec}]: ${stderr}`)
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
