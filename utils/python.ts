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
import type * as path_type from 'path'
import { Readable } from 'stream'

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
export default function python<T>(pathspec: string, args: { kargs?: unknown[], kwargs?: Record<string, unknown> }, callback?: (data: string) => void): Promise<T> {
  if (typeof window !== 'undefined') throw new Error("python is server side only")
  const spawn = typeof window === 'undefined' ? (require('child_process') as typeof child_process_type).spawn : undefined
  const path = typeof window === 'undefined' ? require('path') as typeof path_type : undefined
  return new Promise((resolve, reject) => {
    let stdin: string
    try {
      stdin = JSON.stringify(args)
    } catch (e) {
      throw new ProcessError(`Process input could not be serialized`, null)
    }
    if (typeof spawn === 'undefined') throw new Error("python is server side only")
    if (typeof path === 'undefined') throw new Error("python is server side only")
    const proc = spawn(process.env.PYTHON_BIN || 'python3', [
      path.join(process.env.PYTHON_ROOT || '', 'utils', 'helper.py'),
      pathspec,
    ], { env: { ...process.env } })
    let stdout = ''
    proc.stdout.on('data', (chunk: string) => { stdout += chunk })
    proc.stderr.on('data', callback !== undefined ? callback : (chunk: string) => { console.warn(`[${pathspec}]: ${chunk}`) })
    proc.on('close', (code) => {
      if (code !== 0) {
        reject(new ProcessError(`[${pathspec}]: ${`Process exited with unexpected code ${code}`}`, code))
      } else {
        try {
          resolve(JSON.parse(stdout))
        } catch (e) {
          reject(new ProcessError(`[${pathspec}]: Process output could not be parsed as json`, code))
        }
      }
    })
    proc.stdin.end(stdin)
  })
}


/**
 * @param script: of the form python_file.py:func_in_file
 * @param args: The kwargs to provide to the file
 */
export function pythonStream(pathspec: string, args: { kargs?: unknown[], kwargs?: Record<string, unknown> }): Readable {
  if (typeof window !== 'undefined') throw new Error("python is server side only")
  const spawn = typeof window === 'undefined' ? (require('child_process') as typeof child_process_type).spawn : undefined
  const path = typeof window === 'undefined' ? require('path') as typeof path_type : undefined
  let stdin: string
  try {
    stdin = JSON.stringify(args)
  } catch (e) {
    throw new ProcessError(`Process input could not be serialized`, null)
  }
  if (typeof spawn === 'undefined') throw new Error("python is server side only")
  if (typeof path === 'undefined') throw new Error("python is server side only")
  const proc = spawn(process.env.PYTHON_BIN || 'python3', [
    path.join(process.env.PYTHON_ROOT || '', 'utils', 'helper.py'),
    pathspec,
  ], { env: { ...process.env } })
  let stderr = ''
  proc.stderr.on('data', (chunk: string) => { stderr += chunk })
  proc.on('close', (code) => {
    if (code !== 0) {
      throw new ProcessError(`[${pathspec}]: ${stderr || `Process exited with unexpected code ${code}`}`, code)
    } else {
      if (stderr) {
        console.warn(`[${pathspec}]: ${stderr}`)
      }
    }
  })
  proc.stdin.end(stdin)
  return proc.stdout
}
