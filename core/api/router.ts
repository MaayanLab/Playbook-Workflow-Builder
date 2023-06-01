import { APIRoute } from "@/spec/api";
import getRawBody from 'raw-body'
import * as dict from '@/utils/dict'
import { NotFoundError, UnsupportedMethodError } from "@/spec/error";
import type { NextApiRequest, NextApiResponse } from "next";

function tryJsonParse(v: unknown) {
  try {
    if (typeof v !== 'string') throw new Error()
    return JSON.parse(v)
  } catch (e) {
    return v
  }
}

export default class APIRouter {
  private routes: Record<string, Record<string, APIRoute>> = {}
  private routeExprs: Record<string, RegExp> = {}

  /**
   * Add type-safe API route to in-memory lookup tables
   */
  add(route: APIRoute) {
    if (!(route.path in this.routes)) this.routes[route.path] = {}
    this.routes[route.path][route.method] = route
    if (route.path in this.routeExprs) return
    // convert the path into a RegExp if applicable
    //  /hello/[variable]/world => new RegExp('/hello/(?<variable>[^/]+)/world')
    let expr = /\[(?<arg>.+?)\]/g
    let m = expr.exec(route.path)
    let routeExpr = route.path
    if (!m) return
    while (m !== null) {
      routeExpr = routeExpr.replace(
        m[0],
        m[1].startsWith('..') ? `(?<${m[1].replace(/^\.+/, '')}>.+)` : `(?<${m[1]}>[^/]+)`,
      )
      m = expr.exec(route.path)
    }
    this.routeExprs[route.path] = new RegExp(`^${routeExpr}$`)
  }

  /**
   * Find the API path, resolving path parameter matching when applicable
   */
  find = (path: string) => {
    if (path in this.routes) return { routes: this.routes[path], pathParams: {} }
    const matches: { route: string, m: RegExpExecArray }[] = []
    for (const route in this.routeExprs) {
      const expr = this.routeExprs[route]
      const m = expr.exec(path)
      if (m) matches.push({ route, m })
    }
    if (matches.length > 1) matches.sort((a, b) => a.route.split('/').length - a.route.split('/').length)
    else if (matches.length === 0) return
    const match = matches[0]
    return { routes: this.routes[match.route], pathParams: match.m.groups ?? {} }
  }

  /**
   * Given a NextAPI request, actually call the relevant API and return its response
   */
  route = async (req: NextApiRequest, res: NextApiResponse) => {
    if (!req.url) throw new NotFoundError()
    const url = new URL(req.url, `http://${req.headers.host}`)
    const match = this.find(url.pathname)
    if (typeof match === 'undefined') throw new NotFoundError()
    if (!req.method || !(req.method in match.routes)) throw new UnsupportedMethodError()
    const route = match.routes[req.method]
    const query = route.parameters.parse(
      dict.init(dict.items({ ...match.pathParams, ...Object.fromEntries(url.searchParams.entries()) }).map(({ key, value }) => ({ key, value: tryJsonParse(value) })))
    )
    if (route.method === 'GET') {
      return res.status(200).json(await route.call({ query }, req, res))
    } else if (route.method === 'POST') {
      const rawBody = (await getRawBody(req)).toString()
      const body = route.requestBody.parse(rawBody ? JSON.parse(rawBody) : {})
      return res.status(200).json(await route.call({ query, body }, req, res))
    }
  }
}
