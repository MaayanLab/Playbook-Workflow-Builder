/**
 * Configuration for the browser, read from the environment at runtime
 *  (NEXT_PUBLIC_* variables, e.g. the public url, which is set on deployment)
 */
export function getRuntimeConfig(): Record<string, string> {
  return Object.fromEntries(Object.entries(process.env)
    .filter((entry): entry is [string, string] => entry[0].startsWith('NEXT_PUBLIC_') && entry[1] !== undefined))
}
